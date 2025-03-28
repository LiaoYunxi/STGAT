rm(list=ls())
seed_list=c(10,20,40,60,80,90)

for(seed in seed_list){
  set.seed(seed)
  print(seed)
  source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
  setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/MERFISH/")
  sc.count<-as.matrix(read.table("./sc_data/MERFISH_ID1_0.01_sc_data_transpose.tsv",row.names = 1))
  dim(sc.count)
  st.count<-as.matrix(read.table("./ST_data/MERFISH_ID1_50_data_0.01_transpose.tsv",row.names = 1))
  dim(st.count)
  # saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
  # saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")
  
  lb<-read.table("./sc_data/MERFISH_ID1_sc_label.tsv",sep = "\t",header = T,row.names = 1)
  # rownames(lb) = paste0('X',rownames(lb))
  rownames(lb) = gsub('-','.',rownames(lb))
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
  
  
  
  dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/MERFISH/ID1_0.01_50_DSTG")
  setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/MERFISH/ID1_0.01_50_DSTG")
  
  Convert_Data(count.list,label.list,seed,anova=TRUE)
}

# for(seed in seed_list){
#   set.seed(seed)
#   print(seed)
#   source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/MERFISH/")
#   sc.count<-as.matrix(read.table("./sc_data/MERFISH_ID1_0.01_sc_data_transpose.tsv",row.names = 1))
#   dim(sc.count)
#   st.count<-as.matrix(read.table("./ST_data/MERFISH_ID1_200_data_0.01_transpose.tsv",row.names = 1))
#   dim(st.count)
#   # saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
#   # saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")
#   
#   lb<-read.table("./sc_data/MERFISH_ID1_sc_label.tsv",sep = "\t",header = T,row.names = 1)
#   # rownames(lb) = paste0('X',rownames(lb))
#   rownames(lb) = gsub('-','.',rownames(lb))
#   lb<-as.data.frame(lb[colnames(sc.count),])
#   colnames(lb) =  "subclass"
#   saveRDS(lb,file = "./sc_data/DSTG_sc_label.RDS")
#   
#   intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
#   # id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
#   # intersect.genes =setdiff(intersect.genes ,id)
#   sc.count <- sc.count[intersect.genes,]
#   dim(sc.count)
#   st.count <- st.count[intersect.genes,]
#   dim(st.count)
#   count.list <- list(sc.count,st.count)
#   label.list <- list(data.frame(readRDS("./sc_data/DSTG_sc_label.RDS"),stringsAsFactors=F))
#   
#   
#   
#   dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/MERFISH/ID1_0.01_200_DSTG")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/MERFISH/ID1_0.01_200_DSTG")
#   
#   Convert_Data(count.list,label.list,seed,anova=TRUE)
# }


# st_count=count.list
# st_label=label.list
# anova=TRUE
# data = st_count[[1]]
# label = st_label[[1]]
# nf=2000

# seed=10
# for(seed in seed_list){
#   set.seed(seed)
#   print(seed)
#   source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/seqFISH_plus/Dataset2_seqFISHplus_AllenVIsp/")
#   sc.count<-as.matrix(read.table("./sc_data/sc_data_transpose.tsv",row.names = 1))
#   dim(sc.count)
#   st.count<-as.matrix(read.table("./ST_data/seqFishplus_Cortex_200_data_all_transpose.tsv",row.names = 1))
#   dim(st.count)
#   # saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
#   # saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")
# 
#   # lb<-read.table("./sc_data/sc_label.tsv",sep = "\t",header = T,row.names = 1)
#   # saveRDS(lb,file = "./sc_data/DSTG_sc_label.RDS")
# 
#   intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
#   # id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
#   # intersect.genes =setdiff(intersect.genes ,id)
#   sc.count <- sc.count[intersect.genes,]
#   dim(sc.count)
#   st.count <- st.count[intersect.genes,]
#   dim(st.count)
#   count.list <- list(sc.count,st.count)
#   label.list <- list(data.frame(readRDS("./sc_data/DSTG_sc_label.RDS"),stringsAsFactors=F))
# 
#   dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset2_seqFISHplus_AllenVIsp/DSTG")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset2_seqFISHplus_AllenVIsp/DSTG/")
# 
#   Convert_Data(count.list,label.list,seed,anova=TRUE)
# }
# 
# for(seed in seed_list){
#   set.seed(seed)
#   print(seed)
#   source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/seqFISH_plus/Dataset3_seqFISHplus_sometasensory/")
#   sc.count<-as.matrix(read.table("./sc_data/sc_data_transpose.tsv",row.names = 1))
#   dim(sc.count)
#   st.count<-as.matrix(read.table("./ST_data/seqFishplus_Cortex_200_data_all_transpose.tsv",row.names = 1))
#   dim(st.count)
#   # saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
#   # saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")
#   
#   # lb<-read.table("./sc_data/sc_label.tsv",sep = "\t",header = T,row.names = 1)
#   # lb<-as.data.frame(lb[colnames(sc.count),])
#   # colnames(lb) =  "subclass"
#   # saveRDS(lb,file = "./sc_data/DSTG_sc_label.RDS")
#   
#   intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
#   # id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
#   # intersect.genes =setdiff(intersect.genes ,id)
#   sc.count <- sc.count[intersect.genes,]
#   dim(sc.count)
#   st.count <- st.count[intersect.genes,]
#   dim(st.count)
#   count.list <- list(sc.count,st.count)
#   label.list <- list(data.frame(readRDS("./sc_data/DSTG_sc_label.RDS"),stringsAsFactors=F))
#   
#   dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset3_seqFISHplus_sometasensory/DSTG")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/seqFISH_plus/Dataset3_seqFISHplus_sometasensory/DSTG")
#   
#   Convert_Data(count.list,label.list,seed,anova=TRUE)
# }


# for(seed in seed_list){
#   set.seed(seed)
#   print(seed)
#   source("/data/lyx/software/DSTG-main/DSTG/R_utils.R")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark_data/stereo_seq/zebrafish/")
#   sc.count<-as.matrix(read.table("./sc_data/zebrafish_time18_slice8_sc_data_transpose.tsv",row.names = 1))
#   dim(sc.count)
#   st.count<-as.matrix(read.table("./ST_data/zebrafish_50_data_time18_slice8_transpose.tsv",row.names = 1))
#   dim(st.count)
#   # saveRDS(sc.count,file = "./sc_data/DSTG_sc_data.RDS")
#   # saveRDS(st.count,file = "./ST_data/DSTG_st_data.RDS")
#   
#   # lb<-read.table("./sc_data/zebrafish_time18_slice8_sc_label.tsv",sep = "\t",header = T,row.names = 1)
#   # rownames(lb) = paste0('X',rownames(lb))
#   # rownames(lb) = gsub('-','.',rownames(lb))
#   # lb<-as.data.frame(lb[colnames(sc.count),])
#   # colnames(lb) =  "subclass"
#   # saveRDS(lb,file = "./sc_data/zebrafish_time18_slice8_DSTG_sc_label.RDS")
#   
#   intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
#   # id = c('Foxd2', 'Tlr5', 'Nlrp1c-ps', 'Tifab', 'Gm1966', 'Was', 'Tlr8', 'Slc6a12', 'Cd300ld', 'Tlr13', 'Osm', 'Gsx1', 'Esrp1')
#   # intersect.genes =setdiff(intersect.genes ,id)
#   sc.count <- sc.count[intersect.genes,]
#   dim(sc.count)
#   st.count <- st.count[intersect.genes,]
#   dim(st.count)
#   count.list <- list(sc.count,st.count)
#   label.list <- list(data.frame(readRDS("./sc_data/zebrafish_time18_slice8_DSTG_sc_label.RDS"),stringsAsFactors=F))
#   
#   dir.create("/data/lyx/hubs/SpaTD/stdgcn/benchmark/stereo_seq/zebrafish/zebrafish_time18_slice8_DSTG")
#   setwd("/data/lyx/hubs/SpaTD/stdgcn/benchmark/stereo_seq/zebrafish/zebrafish_time18_slice8_DSTG")
#   
#   Convert_Data(count.list,label.list,seed,anova=TRUE)
# }