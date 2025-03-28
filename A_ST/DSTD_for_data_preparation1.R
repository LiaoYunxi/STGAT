rm(list=ls())
source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
seed=10
set.seed(seed)
###########Dataset1################
# setwd("/data/lyx/software/DSTG-main/DSTG/")
# data<-readRDS("./synthetic_data/example_data.RDS")
# label <-readRDS("./synthetic_data/example_label.RDS")
# d1<-data[[1]]
# l1<-label[[1]]
# l2<-label[[2]]


setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/seqFISH_plus/Dataset1_seqFISHplus_from_stdGCN/")
sc.count<-as.matrix(read.table("./sc_data/sc_data_transpose.tsv",row.names = 1))
dim(sc.count)
st.count<-as.matrix(read.table("./ST_data/ST_data_tanspose.tsv",row.names = 1))
# saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
# saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")

# lb<-read.table("./sc_data/sc_label.tsv",sep = "\t",header = T,row.names = 1)
# saveRDS(lb,file = "./sc_data/DSTG_sc_label.RDS")

intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
sc.count <- sc.count[intersect.genes,]
st.count <- st.count[intersect.genes,]
count.list <- list(sc.count,st.count)
label.list <- list(data.frame(readRDS("./sc_data/DSTG_sc_label.RDS"),stringsAsFactors=F))

dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset1_seqFISHplus_from_stdGCN/DSTG")
setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset1_seqFISHplus_from_stdGCN/DSTG")

Convert_Data(count.list,label.list,seed,anova=TRUE)

# ST_count1<-read.csv("./Infor_Data/ST_count/ST_count_1.csv",row.names = 1)
# ST_count2<-read.csv("./Infor_Data/ST_count/ST_count_2.csv",row.names = 1)
# ST_ref_count1<-read.csv("/data/lyx/software/DSTG-main/DSTG/Infor_Data/ST_count/ST_count_1.csv",row.names = 1)
# # 
# ST_label1<-read.csv("./Infor_Data/ST_label/ST_label_1.csv",row.names = 1)
# ST_ref_label1<-read.csv("/data/lyx/software/DSTG-main/DSTG/Infor_Data/ST_label/ST_label_1.csv",row.names = 1)
# ST_label2<-read.csv("./Infor_Data/ST_label/ST_label_2.csv",row.names = 1)
# ST_ref_label2<-read.csv("/data/lyx/software/DSTG-main/DSTG/Infor_Data/ST_label/ST_label_2.csv",row.names = 1)
# # 
# pseudo_spot <-read.csv("./Datadir/Pseudo_ST1.csv",row.names = 1)
# pseudo_ref_spot <-read.csv("/data/lyx/software/DSTG-main/DSTG/Datadir/Pseudo_ST1.csv",row.names = 1)
# # 
# ST_norm_1 <-read.csv("./Infor_Data/ST_norm/ST_norm_1.csv",row.names = 1)
# ST_ref_norm_1 <-read.csv("/data/lyx/software/DSTG-main/DSTG/Infor_Data/ST_norm/ST_norm_1.csv",row.names = 1)

###########Dataset2################
rm(list=ls())
seed_list=c(10,20,40,60,80,90)
for(seed in seed_list){
  set.seed(seed)
  source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
  setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/seqFISH_plus/Dataset2_seqFISHplus_AllenVIsp/")
  sc.count<-as.matrix(read.table("./sc_data/sc_data_transpose.tsv",row.names = 1))
  dim(sc.count)
  st.count<-as.matrix(read.table("./ST_data/seqFishplus_Cortex_200_data_all_transpose.tsv",row.names = 1))
  dim(st.count)
  # saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
  # saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")
  
  # lb<-read.table("./sc_data/sc_label.tsv",sep = "\t",header = T,row.names = 1)
  # saveRDS(lb,file = "./sc_data/DSTG_sc_label.RDS")
  
  intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
  # id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
  # intersect.genes =setdiff(intersect.genes ,id)
  sc.count <- sc.count[intersect.genes,]
  dim(sc.count)
  st.count <- st.count[intersect.genes,]
  dim(st.count)
  count.list <- list(sc.count,st.count)
  label.list <- list(data.frame(readRDS("./sc_data/DSTG_sc_label.RDS"),stringsAsFactors=F))
  
  dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset2_seqFISHplus_AllenVIsp/DSTG")
  setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset2_seqFISHplus_AllenVIsp/DSTG/")
  
  Convert_Data(count.list,label.list,anova=TRUE)
}

###########Dataset3################
source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/seqFISH_plus/Dataset3_seqFISHplus_sometasensory/")
sc.count<-as.matrix(read.table("./sc_data/sc_data_transpose.tsv",row.names = 1))
dim(sc.count)
st.count<-as.matrix(read.table("./ST_data/seqFishplus_Cortex_200_data_all_transpose.tsv",row.names = 1))
dim(st.count)
saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")

lb<-read.table("./sc_data/sc_label.tsv",sep = "\t",header = T,row.names = 1)
lb<-as.data.frame(lb[colnames(sc.count),])
colnames(lb) =  "subclass"
saveRDS(lb,file = "./sc_data/DSTG_sc_label.RDS")

intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
# id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
# intersect.genes =setdiff(intersect.genes ,id)
sc.count <- sc.count[intersect.genes,]
dim(sc.count)
st.count <- st.count[intersect.genes,]
dim(st.count)
count.list <- list(sc.count,st.count)
label.list <- list(data.frame(readRDS("./sc_data/DSTG_sc_label.RDS"),stringsAsFactors=F))

dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset3_seqFISHplus_sometasensory/DSTG")
setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset3_seqFISHplus_sometasensory/DSTG")

Convert_Data(count.list,label.list,anova=TRUE)
anova=TRUE
st_count=count.list
st_label=label.list
se_obj=tem.t1
clust_vr='subclass'
n=1000

###########zebrafish time12###########
source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/stereo_seq/zebrafish/")
sc.count<-as.matrix(read.table("./sc_data/zebrafish_time12_slice8_sc_data_transpose.tsv",row.names = 1))
dim(sc.count)
st.count<-as.matrix(read.table("./ST_data/zebrafish_50_data_time12_slice8_transpose.tsv",row.names = 1))
dim(st.count)
# saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
# saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")

lb<-read.table("./sc_data/zebrafish_time12_slice8_sc_label.tsv",sep = "\t",header = T,row.names = 1)
rownames(lb) = paste0('X',rownames(lb))
rownames(lb) = gsub('-','.',rownames(lb))
lb<-as.data.frame(lb[colnames(sc.count),])
colnames(lb) =  "subclass"
saveRDS(lb,file = "./sc_data/zebrafish_time12_slice8_DSTG_sc_label.RDS")

intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
# id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
# intersect.genes =setdiff(intersect.genes ,id)
sc.count <- sc.count[intersect.genes,]
dim(sc.count)
st.count <- st.count[intersect.genes,]
dim(st.count)
count.list <- list(sc.count,st.count)
label.list <- list(data.frame(readRDS("./sc_data/zebrafish_time12_slice8_DSTG_sc_label.RDS"),stringsAsFactors=F))

dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/stereo_seq/zebrafish/zebrafish_time12_slice8_DSTG")
setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/stereo_seq/zebrafish/zebrafish_time12_slice8_DSTG")

Convert_Data(count.list,label.list,anova=TRUE)

###########zebrafish time18###########
source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/stereo_seq/zebrafish/")
sc.count<-as.matrix(read.table("./sc_data/zebrafish_time18_slice8_sc_data_transpose.tsv",row.names = 1))
dim(sc.count)
st.count<-as.matrix(read.table("./ST_data/zebrafish_50_data_time18_slice8_transpose.tsv",row.names = 1))
dim(st.count)
# saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
# saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")

lb<-read.table("./sc_data/zebrafish_time18_slice8_sc_label.tsv",sep = "\t",header = T,row.names = 1)
rownames(lb) = paste0('X',rownames(lb))
rownames(lb) = gsub('-','.',rownames(lb))
lb<-as.data.frame(lb[colnames(sc.count),])
colnames(lb) =  "subclass"
saveRDS(lb,file = "./sc_data/zebrafish_time18_slice8_DSTG_sc_label.RDS")

intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
# id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
# intersect.genes =setdiff(intersect.genes ,id)
sc.count <- sc.count[intersect.genes,]
dim(sc.count)
st.count <- st.count[intersect.genes,]
dim(st.count)
count.list <- list(sc.count,st.count)
label.list <- list(data.frame(readRDS("./sc_data/zebrafish_time18_slice8_DSTG_sc_label.RDS"),stringsAsFactors=F))

dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/stereo_seq/zebrafish/zebrafish_time18_slice8_DSTG")
setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/stereo_seq/zebrafish/zebrafish_time18_slice8_DSTG")

Convert_Data(count.list,label.list,anova=TRUE)
