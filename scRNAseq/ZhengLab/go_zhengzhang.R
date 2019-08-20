


setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
#.set_python('/Users/zha8dh/anaconda3/bin/python')

#DATA=readRDS('DATA1.RDS') # Zheng Zhang
#BATCH=readRDS('BATCH1.RDS')


pbmc=readRDS('./pbmc1.RDS')


LABEL=readRDS(file='LABEL.RDS')
pbmc@meta.data$celltype=rep(NA,ncol(pbmc))
pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='NATURE')]=LABEL
#DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=F)
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)
