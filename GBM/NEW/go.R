setwd('/Volumes/Feng/MATT/')



TAG=as.character(read.table('GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.TAG.txt',header=F)[,1])
#a=read.table('GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.txt',row.names=1,header=T,sep='\t',check.names = F)
a=readRDS('GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.RDS')
