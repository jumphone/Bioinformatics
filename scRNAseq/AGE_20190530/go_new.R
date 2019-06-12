
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


#Olmf4Cell

this_gene='Olfm4'
this_index=which(rownames(pbmc)==this_gene)
this_count=pbmc@assays$RNA@counts[this_index,]
Olfm4Cell=which(this_count>0)


this_gene='Ascl2'
this_index=which(rownames(pbmc)==this_gene)
this_count=pbmc@assays$RNA@counts[this_index,]
Ascl2Cell=which(this_count>0)


USED_GENE1=c('Wnt1','Wnt2','Wnt3','Wnt3a','Wnt8a','Wnt8b','Wnt10a','Wnt10b',
             'Ascl2', 'Axin2', 'Olfm4', 'Ctnnb1', 'Ephb2', 'Cd44', 'Myc')

USED_GENE2=c('Wnt4', 'Wnt5a', 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt9a', 'Wnt9b', 'Wnt11')
              


pdf("NEWHEAT.pdf",width=10,height=7)
Olfm4EXP=as.matrix(pbmc@assays$RNA@data[,Olfm4Cell])
Olfm4META=pbmc@meta.data[Olfm4Cell,]
#Olfm4EXP=as.matrix(Olfm4EXP[which(rownames(Olfm4EXP) %in% USED_GENE),])
COL=rep('blue',ncol(Olfm4EXP))
COL[which(Olfm4META$batch=='Aged')]='red'

library('gplots')
heatmap.2(Olfm4EXP[which(rownames(Olfm4EXP) %in% USED_GENE1),],ColSideColors=COL,labCol='',scale=c("none"),dendrogram='none',Rowv=F,Colv=F,trace='none',col=colorRampPalette(c('grey80','red')),margins=c(5,5))

heatmap.2(Olfm4EXP[which(rownames(Olfm4EXP) %in% USED_GENE2),],ColSideColors=COL,labCol='',scale=c("none"),dendrogram='none',Rowv=F,Colv=F,trace='none',col=colorRampPalette(c('grey80','red')),margins=c(5,5))



Ascl2EXP=as.matrix(pbmc@assays$RNA@data[,Ascl2Cell])
Ascl2META=pbmc@meta.data[Ascl2Cell,]
#Olfm4EXP=as.matrix(Olfm4EXP[which(rownames(Olfm4EXP) %in% USED_GENE),])
COL=rep('blue',ncol(Ascl2EXP))
COL[which(Ascl2META$batch=='Aged')]='red'

library('gplots')
heatmap.2(Ascl2EXP[which(rownames(Ascl2EXP) %in% USED_GENE1),],ColSideColors=COL,labCol='',scale=c("none"),dendrogram='none',Rowv=F,Colv=F,trace='none',col=colorRampPalette(c('grey80','red')),margins=c(5,5))

heatmap.2(Ascl2EXP[which(rownames(Ascl2EXP) %in% USED_GENE2),],ColSideColors=COL,labCol='',scale=c("none"),dendrogram='none',Rowv=F,Colv=F,trace='none',col=colorRampPalette(c('grey80','red')),margins=c(5,5))
dev.off()
