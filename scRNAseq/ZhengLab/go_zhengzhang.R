


setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS('./pbmc1.RDS')


LABEL=readRDS(file='LABEL.RDS')
pbmc@meta.data$celltype=rep(NA,ncol(pbmc))
pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='NATURE')]=LABEL
#DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=F)
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)
