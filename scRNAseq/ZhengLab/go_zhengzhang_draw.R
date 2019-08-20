


setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS('./pbmc1.RDS')
pbmc_zhengzhang=readRDS('pbmc_zhengzhang.RDS')

pbmc@meta.data=readRDS('pbmc1_meta.RDS')
pbmc_zhengzhang@meta.data=readRDS('pbmc_zhengzhang_meta.RDS')




DimPlot(pbmc_zhengzhang, reduction.use='umap', group.by='batch', pt.size=0.1,label=T)
