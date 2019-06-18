
library(dplyr)
library(Seurat)
pbmc=readRDS(file='ALL.RDS')

pbmc@meta.data$newall=pbmc@meta.data$all
pbmc@meta.data$newall[which(pbmc@meta.data$all %in% c('TA.Early','TA.G1','TA.G2'))]='TA'

STEM=which( pbmc@meta.data$all %in% c('TA.Early','TA.G1','TA.G2','Enterocyte.Progenitor',
               'Enterocyte.Progenitor.Late') &
           pbmc@assays$RNA@data[which(rownames(pbmc)=='Olfm4'),] >0)
pbmc@meta.data$newall[STEM]='Stem'

FeaturePlot(pbmc,features='Olfm4')


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




STEM=STEM#which( pbmc@meta.data$all %in% c('TA.Early','TA.G1','TA.G2') &  pbmc@assays$RNA@data[which(rownames(pbmc)=='Olfm4'),] >0)
EXP=as.matrix(pbmc@assays$RNA@data[,STEM])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(pbmc@meta.data$batch[STEM]))
OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'PT.cls',sep=' ',quote=F,row.names=F,col.names=F )





USED_GENE1=c('Wnt1','Wnt2','Wnt3','Wnt3a','Wnt8a','Wnt8b','Wnt10a','Wnt10b',
             'Ascl2', 'Axin2', 'Olfm4', 'Ctnnb1', 'Ephb2', 'Cd44', 'Myc')

USED_GENE2=c('Wnt4', 'Wnt5a', 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt9a', 'Wnt9b', 'Wnt11')
          
Olfm4Cell=which(pbmc@meta.data$newall=='Stem')



pdf("NEWNEWHEAT.pdf",width=10,height=7)

Olfm4EXP=as.matrix(pbmc@assays$RNA@data[,Olfm4Cell])
VAR=apply(Olfm4EXP,1,var)
Olfm4EXP=Olfm4EXP[which(VAR>0),]

Olfm4META=pbmc@meta.data[Olfm4Cell,]
COL=rep('blue',ncol(Olfm4EXP))
COL[which(Olfm4META$batch=='Aged')]='red'

library('gplots')
heatmap.2(Olfm4EXP[which(rownames(Olfm4EXP) %in% USED_GENE1),],ColSideColors=COL,labCol='',scale=c("none"),dendrogram='none',Rowv=F,Colv=F,trace='none',col=colorRampPalette(c('grey80','red')),margins=c(5,5))

heatmap.2(Olfm4EXP[which(rownames(Olfm4EXP) %in% USED_GENE2),],ColSideColors=COL,labCol='',scale=c("none"),dendrogram='none',Rowv=F,Colv=F,trace='none',col=colorRampPalette(c('grey80','red')),margins=c(5,10))

dev.off()


