import os
import os.path as osp
import sys
import h5py
import numpy as np
import random
import pandas as pd
import argparse
from tqdm import tqdm
import scipy.sparse
from scipy.sparse import issparse
import scanpy as sc
import anndata as ad

from process_img import *
from utils import *

import torch
import torch.nn as nn
import torch.nn.functional as F
from PIL import Image
from torch.utils.data import Dataset, DataLoader, random_split
from torchvision import transforms

sys.path.append("/data/lyx/hubs/Cpath/src/")
from models import *
import time
from PIL import Image

bin_size = 100
crop_size = 150
output='/data/lyx/hubs/Cpath/hest_data/image_out/'
name = "TENX94"
name = "".join(name)
label_txt_file=str(name+'.h5')
adata_file=str(name+'.h5ad')
tile_path = os.path.join(output,name)
current_path = "/data/lyx/software/StereoMMv1/"
data_path = "/data/lyx/hubs/Cpath/hest_data/"
Image.MAX_IMAGE_PIXELS = None

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu") # device object

with h5py.File(osp.join(data_path,label_txt_file), "r") as f:
    # 获取数据集对象
    img_dataset = f["img"]
    # 将数据集读取为NumPy数组
    img_data = img_dataset[:]  # 或 np.array(img_dataset)
    barcode_set = f["barcode"]
    # 将数据集读取为NumPy数组
    barcode = barcode_set[:]  # 或 np.array(img_dataset)

string_data = [byte_data.decode('utf-8') for byte_data in barcode.flatten().tolist()]
    
print("数据形状:", img_data.shape)  # 应为 (6078, 224, 224, 3)
print("数据类型:", img_data.dtype)  # 应为 uint8

model_com = models.resnet50(pretrained=False)
pth_path = os.path.join(current_path, "torch_pths/resnet50-19c8e357.pth")
model_com.load_state_dict(torch.load(pth_path))
num_features = model_com.fc.in_features

img_model = torch.nn.Sequential(*list(model_com.children())[:-1])

img_data = torch.from_numpy(img_data).permute(0, 3, 1, 2).float()  # 转换为浮点型

train_transform = transforms.Compose([
    transforms.RandomHorizontalFlip(),
    transforms.RandomResizedCrop(224, scale=(0.8, 1.0)),
    transforms.ColorJitter(brightness=0.2, contrast=0.2),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])

tgthighres = 'a5'

ckpt_path = '/data/lyx/software/scFoundation/model/models/models.ckpt'
pretrainmodel,pretrainconfig = load_model_frommmf(ckpt_path,"cpu","rde")

adata = ad.read_h5ad(f"./hest_data/{name}.h5ad")
adata =adata[string_data,:].copy()
adata.write_h5ad(f"./hest_data/{name}.h5ad")

# 初始化存储列表
x_list = []
# 处理每个细胞
for i in tqdm(range(adata.shape[0])):
    # 获取原始counts数据（假设X矩阵存储原始counts）
    row = adata.X[i, :]
    
    # 处理稀疏矩阵
    if not isinstance(row, np.ndarray):
        row = row.toarray().flatten()
    else:
        row = row.flatten()
    
    # 计算归一化后的log1p值
    totalcount = row.sum()
    if totalcount == 0:
        totalcount = 1e-5  # 避免除零错误
    normalized = (row / totalcount) * 1e4
    tmpdata = np.log1p(normalized).tolist()
    
    # 构建特征张量
    log_total = np.log10(totalcount)
    pretrain_gene_x = torch.tensor(
        tmpdata + [log_total + float(tgthighres[1:]), log_total]
    ).unsqueeze(0)
    
    
    # x, position_gene_ids,x_padding, encoder_data_labels=getEncoderData(pretrain_gene_x,pretrainconfig)
    # x_list.append((x,position_gene_ids,img_emb_all[i,:]))#x_padding
    
    # 生成基因位置ID
    data_gene_ids = torch.arange(19266, device=pretrain_gene_x.device).repeat(pretrain_gene_x.shape[0], 1)
    
    # 数据打包
    value_labels = pretrain_gene_x > 0
    x, x_padding = gatherData(pretrain_gene_x, value_labels, pretrainconfig['pad_token_id'])
    position_gene_ids, _ = gatherData(data_gene_ids, value_labels, pretrainconfig['pad_token_id'])
    
    # 存储结果
    x_list.append((x,position_gene_ids,img_data[i,:,:,:]))


random.seed(0)
np.random.seed(0)  # numpy random generator

torch.manual_seed(0)
torch.cuda.manual_seed_all(0)

torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

model = ContrastModel(768, img_model, pretrainmodel,pretrainconfig,frozenmore=False)

source_state_dict = pretrainmodel.state_dict()
target_state_dict = model.state_dict()

encoder_keys = [k for k in source_state_dict.keys() if k.startswith("encoder.")]
encoder_keys_target = [k for k in target_state_dict.keys() if k.startswith("encoder.")]
encoder_keys = [i for i in encoder_keys if i in encoder_keys_target]

source_state_dict = {k: v for k, v in source_state_dict.items() if k in encoder_keys}

for k,v in source_state_dict.items():
    target_state_dict[k]=v

model.load_state_dict(target_state_dict, strict=False)
model.build()

model = model.to(device)

dataset = CustomDataset(x_list,transform=train_transform )
dataloader = DataLoader(
    dataset,
    batch_size=32,
    collate_fn=collate_fn,
    shuffle=True
)

epochs = 100
warmup_steps = 1000
clip_grad = 2.0

logit_scale = nn.Parameter(torch.tensor(np.log(1/0.07)))
logit_scale.data = torch.clamp(logit_scale, min=np.log(1e-4), max=np.log(100)).detach()

# 初始化优化器
optimizer = torch.optim.AdamW(
    get_parameter_groups(model,logit_scale),
    betas=(0.9, 0.98),  # 调整beta参数
    eps=1e-6
)

# 添加学习率调度器
# 组合调度策略
scheduler = torch.optim.lr_scheduler.SequentialLR(
    optimizer,
    schedulers=[
        LinearLR(optimizer, start_factor=0.01, total_iters=1000),  # Warmup
        CosineAnnealingLR(optimizer, T_max=epochs*len(dataloader))  # 余弦退火
    ],
    milestones=[1000]
)

scaler = torch.cuda.amp.GradScaler()


for epoch in range(epochs):
    model.train()
    total_loss = 0
    
    for step, batch in enumerate(dataloader):
        optimizer.zero_grad()
        
        # 混合精度前向
        with torch.autocast(device_type='cuda', dtype=torch.float16):
            tx_emb, imge_emb = model(batch['x'].to(device), batch['x_padding'].to(device),
                                     batch['position_gene_ids'].to(device), batch['img'].to(device))
            logits = (imge_emb @ tx_emb.T) * logit_scale.exp()
            loss = contrastive_loss(logits)
        
        # 反向传播
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        
        # 梯度裁剪与更新
        torch.nn.utils.clip_grad_norm_(model.parameters(), clip_grad)
        scaler.step(optimizer)
        scaler.update()
        scheduler.step()
        
        # 记录损失
        total_loss += loss.item()
        
        # 每100步打印日志
        if step % 100 == 0:
            lr = optimizer.param_groups[0]['lr']
            print(f"Epoch {epoch} Step {step} | Loss: {loss.item():.4f} | LR: {lr:.2e}")

    # 保存检查点
    if epoch % 10 == 0:
        torch.save({
            'model': model.state_dict(),
            'optimizer': optimizer.state_dict(),
            'epoch': epoch
        }, f"checkpoint_epoch{epoch}.pth")

# 最终模型保存
torch.save(model.state_dict(), "final_model.pth")