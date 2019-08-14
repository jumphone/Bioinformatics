

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('/Volumes/Feng/MATT/')






library(Seurat)
library(dplyr)




D1=read.delim('matt_combined.txt',header=T,row.names=1,sep='\t')
D2=read.delim('matt_combined2.txt',header=T,row.names=1,sep='\t')




pbmc1 <- CreateSeuratObject(counts =D1, project = "D1", min.cells = 0, min.features = 0)
pbmc1 <- NormalizeData(object = pbmc1, normalization.method = "LogNormalize",  scale.factor = 10000)

pbmc2 <- CreateSeuratObject(counts =D2, project = "D2", min.cells = 0, min.features = 0)
pbmc2 <- NormalizeData(object = pbmc2, normalization.method = "LogNormalize",  scale.factor = 10000)


pbmc=readRDS(file='newpbmc.RDS')

VF=VariableFeatures(object = pbmc)


ND1=as.matrix(pbmc1@assays$RNA@data[which(rownames(pbmc1) %in% VF),])
ND2=as.matrix(pbmc2@assays$RNA@data[which(rownames(pbmc2) %in% VF),])

ND=as.matrix(pbmc@assays$RNA@data[which(rownames(pbmc) %in% VF),])

PC1=pbmc@reductions$pca@cell.embeddings[,1]

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


exp_sc_mat=ND2
exp_ref_mat=ND

out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
tag=.get_tag_max(out)


