import sys
import random,os
import numpy as np
import pandas as pd
import argparse
import torch
import time
from PIL import Image
from tqdm import tqdm
import scipy.sparse
from scipy.sparse import issparse
import scanpy as sc
import anndata as ad
# from select_model import *

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset,DataLoader
from torch.nn.utils.rnn import pad_sequence
from torch.optim.lr_scheduler import CosineAnnealingLR, LinearLR
from torchvision import transforms
from select_model import select_model

def gatherData(data, labels, pad_token_id):
    value_nums = labels.sum(1)
    max_num = max(value_nums)


    fake_data = torch.full((data.shape[0], max_num), pad_token_id,
                           device=data.device)
    data = torch.hstack([data, fake_data])

    fake_label = torch.full((labels.shape[0], max_num), 1,
                            device=labels.device)
    none_labels = ~labels
    labels = labels.float()
    labels[none_labels] = torch.tensor(-float('Inf'), device=labels.device)

    tmp_data = torch.tensor([(i + 1) * 20000 for i in range(labels.shape[1], 0, -1)], device=labels.device)
    labels += tmp_data

    labels = torch.hstack([labels, fake_label])

    fake_label_gene_idx = labels.topk(max_num).indices

    new_data = torch.gather(data, 1, fake_label_gene_idx)

    padding_labels = (new_data == pad_token_id)

    return new_data, padding_labels

def getEncoderData(data, config):
    encoder_data_labels = data > 0
    encoder_data, encoder_data_padding = gatherData(data, encoder_data_labels,
                                                    config['pad_token_id'])
    data_gene_ids = torch.arange(data.shape[1], device=data.device).repeat(data.shape[0], 1)
    encoder_position_gene_ids, _ = gatherData(data_gene_ids, encoder_data_labels,
                                                config['pad_token_id'])
    encoder_position_gene_ids[encoder_data_padding] = config["seq_len"]
    
    return encoder_data, encoder_position_gene_ids, encoder_data_padding, encoder_data_labels

def convertconfig(ckpt):
    newconfig = {}
    newconfig['config']={}
    model_type = ckpt['config']['model']
    
    for key, val in ckpt['config']['model_config'][model_type].items():
        newconfig['config'][key]=val
        
    for key, val in ckpt['config']['dataset_config']['rnaseq'].items():
        newconfig['config'][key]=val
        
    if model_type == 'performergau_resolution':
        model_type = 'performer_gau'
    
    import collections
    d = collections.OrderedDict()
    for key, val in ckpt['state_dict'].items():
        d[str(key).split('model.')[1]]=val
        
    newconfig['config']['model_type']=model_type
    newconfig['model_state_dict']=d
    newconfig['config']['pos_embed']=False
    newconfig['config']['device']='cuda'
    return newconfig

def load_model_frommmf(best_ckpt_path, gpu="cpu",key='gene'):
    model_data = torch.load(best_ckpt_path,map_location=gpu)
    model_data = model_data[key]
    model_data = convertconfig(model_data)

    if not model_data.__contains__('config'):
        print('***** No config *****')
        config={}
        config['model_type']='flash_all'
    else:
        config=model_data['config']
        print(config)
    if not config.__contains__('qv_dim'):
        if config['model'] != 'mae_autobin':
            if config.__contains__('dim_head'):
                config['qv_dim']=config['dim_head']
            else:
                print('***** No qv_dim ***** set 64')
                config['qv_dim']= 64
    if not config.__contains__('ppi_edge'):
        config['ppi_edge']=None
    model = select_model(config)
    model_state_dict = model_data['model_state_dict']    
    model.load_state_dict(model_state_dict)
    return model,config

def get_parameter_groups(model,logit_scale,pretrain_lr=1e-5, base_lr=1e-4, img_lr=5e-5):
    params = []
    
    # 预训练模型参数组（更低学习率）
    params.append({
        'params': [p for p in model.encoder.parameters() if p.requires_grad],
        'lr': pretrain_lr,
        'weight_decay': 0.0  # 通常预训练模型不适用权重衰减
    })
    
    # 文本特征融合层（中等学习率）
    text_params = list(model.fc1.parameters()) + list(model.fc2.parameters())
    params.append({
        'params': text_params,
        'lr': base_lr,
        'weight_decay': 0.05
    })
    
    # 图像处理层（独立学习率）
    img_params = list(model.img_fc1.parameters()) + \
                list(model.img_fc2.parameters()) 
    params.append({
        'params': img_params,
        'lr': img_lr,
        'weight_decay': 0.05
    })
    
    # 温度参数（独立配置）
    params.append({
        'params': [logit_scale],
        'lr': 0.001,
        'weight_decay': 0.0
    })
    
    return params

def contrastive_loss(logits, smooth=0.1):
    batch_size = logits.size(0)
    labels = torch.arange(batch_size, device=logits.device)
    loss_i = F.cross_entropy(logits, labels, label_smoothing=smooth)
    loss_t = F.cross_entropy(logits.T, labels, label_smoothing=smooth)
    return (loss_i + loss_t) / 2


class CustomDataset(Dataset):
    def __init__(self, data_list,transform = None):
        self.data = data_list  # 假设每个元素是包含x, pos, pad, label的元组
        self.transform = transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        x,pos,image = self.data[index]
        
        if self.transform:
            image = self.transform(image)
        # 转换为张量
        return (
            torch.tensor(x),
            torch.tensor(pos),
            image
        )
    
def collate_fn(batch):
    xs, position_gene_ids_list, pads_list, img = zip(*batch)
    
    # Calculate max_len based on both x and pos lengths
    max_len_x = max(len(x) for x in xs)
    max_len_pos = max(len(pos) for pos in position_gene_ids_list)
    max_len = max(max_len_x, max_len_pos)  # Use the larger of the two
    
    # Pad input sequences (x)
    padded_x = torch.zeros(len(xs), max_len, dtype=torch.long)
    x_padding = torch.zeros(len(xs), max_len, dtype=torch.bool)
    for i, x in enumerate(xs):
        padded_x[i, :len(x)] = x
        x_padding[i, :len(x)] = pads_list[i]
    
    # Pad positional IDs (pos)
    padded_position = torch.zeros(len(xs), max_len, dtype=torch.long)
    for i, pos_ids in enumerate(position_gene_ids_list):
        padded_position[i, :len(pos_ids)] = pos_ids
    
    # Process images
    img = torch.stack(img, dim=0)
    
    return {
        'x': padded_x,
        'x_padding': x_padding,
        'position_gene_ids': padded_position,
        'img': img
    }

class ContrastModel(nn.Module):
    def __init__(self, input_dim, img_model, pretrainmodel,pretrainconfig, frozenmore=True):
        super().__init__()
        self.frozenmore = frozenmore
        # 文本编码部分
        self.token_emb = pretrainmodel.token_emb
        self.pos_emb = pretrainmodel.pos_emb
        # self.encoder = pretrainmodel.encoder
        
        self.encoder = pytorchTransformerModule(
            max_seq_len = 15000,
            dim = 768,
            depth=12,
            heads=12,
            ff_mult=4,
            norm_first=False
        )
        
        # 文本特征融合层
        self.fc1 = nn.Sequential(
            nn.Linear(input_dim*2, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU()
        )
        self.fc2 = nn.Sequential(
            nn.Linear(input_dim*2, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU()
        )
        
        # 图像编码部分
        self.img_model = img_model

        self.img_fc1 = nn.Sequential(
            # nn.Linear(2048, 1024),  # 假设图像特征维度为2048
            nn.AdaptiveAvgPool1d(1024),
            nn.BatchNorm1d(1024),
            nn.ReLU()
        )
        self.img_fc2 = nn.Sequential(
            nn.Linear(2048, 1024),
            nn.BatchNorm1d(1024),
            nn.ReLU()
        )

        # self.img_proj = nn.Linear(2048, 1024)  # 投影到对比空间
    def build(self):
        if self.frozenmore:
            for _,p in self.token_emb.named_parameters():
                p.requires_grad = False
            for _,p in self.pos_emb.named_parameters():
                p.requires_grad = False
            print('self.pos_emb and self.token_emb also frozen')

            for param in self.img_model.parameters():
                param.requires_grad = False
        
        for na, param in self.encoder.named_parameters():
            param.requires_grad = False
        for na, param in self.encoder.transformer_encoder[-2].named_parameters():
            print('self.encoder.transformer_encoder ',na,' have grad')
            param.requires_grad = True
            
    def forward(self, x,x_padding,position_gene_ids, img):
        x_img = self.img_model(img)
        x = self.token_emb(torch.unsqueeze(x, 2).float(), output_weight=0)

        position_emb = self.pos_emb(position_gene_ids)

        x += position_emb
        geneemb = self.encoder(x, x_padding)
        # print("encoder:",torch.isnan(geneemb).any())
        
        # 多层级特征融合
        geneemb1 = geneemb[:, -1, :]
        geneemb2 = geneemb[:, -2, :]
        geneemb3, _ = torch.max(geneemb[:, :-2, :], dim=1)
        geneemb4 = torch.mean(geneemb[:, :-2, :], dim=1)
        
        x1 = torch.cat([geneemb1, geneemb2], dim=1)
        x2 = torch.cat([geneemb3, geneemb4], dim=1)
        tx_feat = torch.cat([self.fc1(x1), self.fc2(x2)], dim=1)
        
        # 图像编码
        # self.img_model.eval()
        x_img = self.img_model(img)
        x_img = torch.cat([self.img_fc1(x_img),self.img_fc2(x_img)], dim=1)
        # img_feat = self.img_proj(img_feat)
        
        return F.normalize(tx_feat, p=2, dim=1), F.normalize(x_img, p=2, dim=1)


# class ContrastModel(nn.Module):
#     def __init__(self, input_dim, img_model, pretrainmodel,pretrainconfig, frozenmore=True):
#         super().__init__()
#         self.frozenmore = frozenmore
#         # 文本编码部分
#         self.token_emb = pretrainmodel.token_emb
#         self.pos_emb = pretrainmodel.pos_emb
#         # self.encoder = pretrainmodel.encoder
        
#         self.encoder = pytorchTransformerModule(
#             max_seq_len = 15000,
#             dim = 768,
#             depth=12,
#             heads=12,
#             ff_mult=4,
#             norm_first=False
#         )
        
#         # 文本特征融合层
#         self.fc1 = nn.Sequential(
#             nn.Linear(input_dim*2, 1024),
#             nn.BatchNorm1d(1024),
#             nn.ReLU()
#         )
#         self.fc2 = nn.Sequential(
#             nn.Linear(input_dim*2, 1024),
#             nn.BatchNorm1d(1024),
#             nn.ReLU()
#         )
        
#         # 图像编码部分
#         self.img_model = img_model

#         self.img_fc1 = nn.Sequential(
#             # nn.Linear(2048, 1024),  # 假设图像特征维度为2048
#             nn.AdaptiveAvgPool1d(1024),
#             nn.BatchNorm1d(1024),
#             nn.ReLU()
#         )
#         self.img_fc2 = nn.Sequential(
#             nn.Linear(2048, 1024),
#             nn.BatchNorm1d(1024),
#             nn.ReLU()
#         )

#         # self.img_proj = nn.Linear(2048, 1024)  # 投影到对比空间
#     def build(self):
#         if self.frozenmore:
#             for _,p in self.token_emb.named_parameters():
#                 p.requires_grad = False
#             for _,p in self.pos_emb.named_parameters():
#                 p.requires_grad = False
#             print('self.pos_emb and self.token_emb also frozen')
        
#         for na, param in self.encoder.named_parameters():
#             param.requires_grad = False
#         for na, param in self.encoder.transformer_encoder[-2].named_parameters():
#             print('self.encoder.transformer_encoder ',na,' have grad')
#             param.requires_grad = True
        
#         for param in self.img_model.parameters():
#             param.requires_grad = False
            
#     def forward(self, x,x_padding,position_gene_ids, img_feat):
#         x = self.token_emb(torch.unsqueeze(x, 2).float(), output_weight=0)

#         position_emb = self.pos_emb(position_gene_ids)

#         x += position_emb
#         geneemb = self.encoder(x, x_padding)
#         # print("encoder:",torch.isnan(geneemb).any())
        
#         # 多层级特征融合
#         geneemb1 = geneemb[:, -1, :]
#         geneemb2 = geneemb[:, -2, :]
#         geneemb3, _ = torch.max(geneemb[:, :-2, :], dim=1)
#         geneemb4 = torch.mean(geneemb[:, :-2, :], dim=1)
        
#         x1 = torch.cat([geneemb1, geneemb2], dim=1)
#         x2 = torch.cat([geneemb3, geneemb4], dim=1)
#         tx_feat = torch.cat([self.fc1(x1), self.fc2(x2)], dim=1)
        
#         # 图像编码
#         self.img_model.eval()
#         # img_feat = self.img_model(img_emb)
#         img_feat = torch.cat([self.img_fc1(img_feat),self.img_fc2(img_feat)], dim=1)
#         # img_feat = self.img_proj(img_feat)
        
#         return F.normalize(tx_feat, p=2, dim=1), F.normalize(img_feat, p=2, dim=1)

class pytorchTransformerModule(nn.Module):
    def __init__(self,
                 max_seq_len,
                 dim,
                 depth,
                 heads,
                 ff_mult=4,
                 norm_first=False,
                 ):
        super(pytorchTransformerModule, self).__init__()

        self.max_seq_len = max_seq_len
        self.depth = depth
        layers = []
        for i in range(depth):
            layers.append(nn.TransformerEncoderLayer(d_model=dim, nhead=heads,
                                                     dim_feedforward=dim * ff_mult,
                                                     batch_first=True,
                                                     norm_first=norm_first,
                                                     #activation="gelu",
                                                     ))

        self.transformer_encoder = nn.ModuleList(layers)
        self.norm = nn.LayerNorm(dim)

    def forward(self, x, padding_mask):
        b, n, _, device = *x.shape, x.device
        assert n <= self.max_seq_len, f'sequence length {n} must be less than the max sequence length {self.max_seq_len}'

        # x get encodings [B, N, D] , batch_first is True
        count = 0
        for mod in self.transformer_encoder:
            x = mod(x, src_key_padding_mask=padding_mask) # , src_mask=mask, src_key_padding_mask=src_key_padding_mask)
            # print(f"layer {str(count)}:",torch.isnan(x).any())
        # x = self.transformer_encoder(x)
        x = self.norm(x)

        return x

# class ImageEncoder(nn.Module):
#     def __init__(self,img_model):
#         super().__init__()
#         self.model = img_model
#         self.pooling = nn.AdaptiveAvgPool1d(1024)
#         self.fc = nn.Linear(in_features=2048, out_features=1024, bias=True)
#     def forward(self, x):
#         x = self.model(x)
#         x = torch.cat(self.pooling(x),self.fc(x),axis=1)
#         return F.normalize(x, p=2, dim=1)

# class TranscriptomeEncoder(nn.Module):
#     def __init__(self, input_dim,pretrainmodel):
#         super().__init__()
#         self.token_emb = pretrainmodel.token_emb
#         self.pos_emb = pretrainmodel.pos_emb
#         self.encoder = pretrainmodel.encoder
#         self.fc1 = nn.Linear(in_features=input_dim*2, out_features=1024, bias=True)
#         self.fc2 = nn.Linear(in_features=input_dim*2, out_features=1024, bias=True)  
#     def forward(self, x,position_gene_ids,x_padding):
#         x = self.token_emb(x, output_weight = 0)
#         position_emb = self.pos_emb(position_gene_ids)
#         x += position_emb
#         geneemb = self.encoder(x,x_padding)
#         geneemb1 = geneemb[:,-1,:]
#         geneemb2 = geneemb[:,-2,:]
#         geneemb3, _ = torch.max(geneemb[:,:-2,:], dim=1)
#         geneemb4 = torch.mean(geneemb[:,:-2,:], dim=1)
#         x1 = torch.concat([geneemb1,geneemb2],axis=1)
#         x2 = torch.concat([geneemb3,geneemb4],axis=1)
#         x = torch.concat([self.fc1(x1),self.fc1(x2)],axis=1)
#         return F.normalize(x, p=2, dim=1)