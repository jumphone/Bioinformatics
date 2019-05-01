source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')


library(Seurat)
library(dplyr)
library(Matrix)
load('HDN_Wkspace.RData')

pbmc=tumor2
pbmc.raw.data=as.matrix(pbmc@raw.data[,which(colnames(pbmc@raw.data) %in% colnames(pbmc@scale.data))])
pbmc.data=as.matrix(pbmc@scale.data)

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
ONE=.data2one(pbmc.raw.data, rownames(pbmc.data), CPU=4, PCNUM=50, SEED=123,  PP=30)
#ONE=readRDS('ONE.RDS')


OUT=getBIN(ONE)
BIN=OUT$BIN
BINTAG=OUT$TAG

pbmc@meta.data$bin=BINTAG
png('ID.png',width=1200,height=1000)
DimPlot(pbmc,group.by='bin',reduction.use='umap',do.label=T)
dev.off()



EXP=pbmc.data
LR=read.table('RL_mouse.txt',header=T,sep='\t')

MEAN=getMEAN(EXP, LR)

PMAT=getPMAT(EXP, LR, BIN, MEAN)
CMAT=getCMAT(EXP,LR,PMAT)

library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))


OUT=getPAIR(CMAT)

PAIR=OUT$PAIR
SCORE=OUT$SCORE
RANK=OUT$RANK

VEC=pbmc@dr$umap@cell.embeddings
CPlot(VEC,PAIR[1:10,],BINTAG)


ORITAG=as.character(pbmc@ident)
ORITAG[which(ORITAG=='Pericyte/\nFibroblast')]='Pericyte Fibroblast'
NET=getNET(PAIR, BINTAG,ORITAG )



CN=getCN(NET)
DP=DPlot(NET, CN,COL=3)


par(mfrow=c(2,2))
LP=LPlot('Pericyte Fibroblast','Tumor Cells',NET,PMAT,SEED=12345)
LP=LPlot('Pericyte Fibroblast','Proliferating Cells',NET,PMAT,SEED=12345)
LP=LPlot('Oligodendrocytes','Tumor Cells',NET,PMAT,SEED=12345)
LP=LPlot('Oligodendrocytes','Proliferating Cells',NET,PMAT,SEED=12345)
