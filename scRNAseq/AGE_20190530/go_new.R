
library(dplyr)
library(Seurat)
pbmc=readRDS(file='ALL.RDS')

#DimPlot(pbmc, group.by='all',label=T)
pbmc@meta.data$young=pbmc@meta.data$all
pbmc@meta.data$young[which(pbmc@meta.data$batch=='Aged')]=NA
pbmc@meta.data$aged=pbmc@meta.data$all
pbmc@meta.data$aged[which(pbmc@meta.data$batch=='Young')]=NA

AC=which(pbmc@meta.data$batch=='Aged')
YC=which(pbmc@meta.data$batch=='Young')


pdf("YA.pdf",width=10,height=7)
DimPlot(pbmc, cells=colnames(pbmc)[AC], group.by='all',label=T)
DimPlot(pbmc, cells=colnames(pbmc)[YC], group.by='all',label=T)
dev.off()
