
library(dplyr)
library(Seurat)
pbmc=readRDS(file='ALL.RDS')

pbmc@meta.data$newall=pbmc@meta.data$all
pbmc@meta.data$newall[which(pbmc@meta.data$all %in% c('TA.Early','TA.G1','TA.G2'))]='TA'

STEM=which( pbmc@meta.data$newall %in% c('TA') &  pbmc@assays$RNA@data[which(rownames(pbmc)=='Olfm4'),] >0)
pbmc@meta.data$newall[STEM]='Stem'




#DimPlot(pbmc, group.by='all',label=T)
pbmc@meta.data$young=pbmc@meta.data$newall
pbmc@meta.data$aged=pbmc@meta.data$newall


pbmc@meta.data$young[which(pbmc@meta.data$batch=='Aged')]=NA
pbmc@meta.data$aged[which(pbmc@meta.data$batch=='Young')]=NA

AC=which(pbmc@meta.data$batch=='Aged')
YC=which(pbmc@meta.data$batch=='Young')


pdf("YA_STEM.pdf",width=10,height=7)
DimPlot(pbmc, cells=colnames(pbmc), group.by='newall',label=T)
DimPlot(pbmc, cells=colnames(pbmc)[AC], group.by='newall',label=T)
DimPlot(pbmc, cells=colnames(pbmc)[YC], group.by='newall',label=T)
dev.off()



OUT=t(table(pbmc@meta.data$batch,pbmc@meta.data$newall))
write.table(OUT,file='BSTAT.txt',sep='\t',row.names=T,col.names=T,quote=F)
