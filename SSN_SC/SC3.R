library(SingleCellExperiment)
library(SC3)
library(scater)
exp_data=read.table('test0.2.data.result.b',header=T,row.names=1)
exp_data[is.na(exp_data)]=0
SUM_CELL=apply(abs(exp_data),2,sum)
SUM_GENE=apply(abs(exp_data),1,sum)

exp_data=exp_data[which(SUM_GENE>0),which(SUM_CELL>0)]

sce <- SingleCellExperiment(assays = list(counts=as.matrix(exp_data),logcounts=as.matrix(exp_data)), colData=colnames(exp_data))
plotPCA(sce)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3(sce, ks = 2:4, biology = TRUE)

