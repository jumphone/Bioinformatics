source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')


library(Seurat)
library(dplyr)
library(Matrix)


load('Seurat_EXP_cluster.Robj')

pbmc.raw.data=as.matrix(EXP_cluster@raw.data[,which(colnames(EXP_cluster@raw.data) %in% colnames(EXP_cluster@scale.data))])
pbmc.data=as.matrix(EXP_cluster@scale.data)

#source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
#ONE=.data2one(pbmc.raw.data, rownames(pbmc.data), CPU=4, PCNUM=50, SEED=123,  PP=30)
ONE=readRDS('ONE.RDS')


BIN=readRDS('BIN.RDS')


BINTAG=rep(NA,length(ONE))
i=1
while(i<=ncol(BIN)){
BINTAG[BIN[,i]]=i
i=i+1
}


pbmc=EXP_cluster
pbmc@meta.data$bin=BINTAG
png('ID.png',width=1200,height=1000)
DimPlot(pbmc,group.by='bin',reduction.use='tsne',do.label=T)
dev.off()



EXP=pbmc.data
LR=read.table('RL_mouse.txt',header=T,sep='\t')

MEAN=readRDS('MEAN.RDS')
PMAT=readRDS('PMAT.RDS')
CMAT=readRDS('CMAT.RDS')

library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))


OUT=getPAIR(CMAT)

PAIR=OUT$PAIR
SCORE=OUT$SCORE
RANK=OUT$RANK

VEC=pbmc@dr$tsne@cell.embeddings
CPlot(VEC,PAIR[1:200,],BINTAG)







ORITAG=rep('NA',length(ONE))
ORITAG[which(pbmc@meta.data$orig.ident=='S1.5MD' & BINTAG %in% c(73,74))]='1.5mo_73,74_Fibroblasts'
ORITAG[which(pbmc@meta.data$orig.ident=='S1.5MD' & BINTAG %in% c(73,74))]='1.5mo_30,31,32,33_Tcell'

as.character(pbmc@ident)




ORITAG[which(ORITAG=='Pericyte/\nFibroblast')]='Pericyte Fibroblast'
NET=getNET(PAIR[1:200], BINTAG,ORITAG )



CN=getCN(NET)
DP=DPlot(NET, CN,COL=3)


par(mfrow=c(2,2))
LP=LPlot('Pericyte Fibroblast','Tumor Cells',NET,PMAT,SEED=12345)
LP=LPlot('Pericyte Fibroblast','Proliferating Cells',NET,PMAT,SEED=12345)
LP=LPlot('Oligodendrocytes','Tumor Cells',NET,PMAT,SEED=12345)
LP=LPlot('Oligodendrocytes','Proliferating Cells',NET,PMAT,SEED=12345)
