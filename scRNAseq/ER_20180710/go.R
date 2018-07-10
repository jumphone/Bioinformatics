
library(Seurat)
library(dplyr)
library(Matrix)

PCNUM=40
PCUSE=1:35
RES=2
raw_data=read.table('EditingMatrix.tsv',header=T,row.names=NULL,check.names=F)
raw_data=t(raw_data)
head(rownames(raw_data),n=20)

dim(raw_data)
exp_data=raw_data[c(11:length(raw_data[,1])),]
exp_data=as.matrix(exp_data)
exp_data=apply(exp_data, 2, as.numeric)
colnames(exp_data)=raw_data[1,]
rownames(exp_data)=rownames(raw_data)[c(11:length(raw_data[,1]))]

EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes= 0)

head(rownames(raw_data),n=20)
EXP <- AddMetaData(object = EXP, metadata = raw_data[2,], col.name = "sample_type")
EXP <- AddMetaData(object = EXP, metadata = raw_data[3,], col.name = "tumor_stage")
EXP <- AddMetaData(object = EXP, metadata = raw_data[5,], col.name = "gender")
EXP <- AddMetaData(object = EXP, metadata = raw_data[6,], col.name = "race")
EXP <- AddMetaData(object = EXP, metadata = raw_data[7,], col.name = "ethnicity")
EXP <- AddMetaData(object = EXP, metadata = as.numeric(raw_data[8,]), col.name = "days_to_birth")
EXP <- AddMetaData(object = EXP, metadata = as.numeric(raw_data[9,]), col.name = "days_to_death")
EXP <- AddMetaData(object = EXP, metadata = as.numeric(raw_data[10,]), col.name = "project_id")




