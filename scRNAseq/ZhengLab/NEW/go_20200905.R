setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
library(Seurat)
pbmc=readRDS('../pbmc.final.RDS')

wt.pbmc=readRDS('new.wt.pbmc.rds')
ko.pbmc=readRDS('new.ko.pbmc.rds')

HQRef=readRDS('HQRef.rds')

