source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')


library(Seurat)
library(dplyr)
library(Matrix)
load('HDN_Wkspace.RData')

pbmc=tumor2
pbmc.raw.data=as.matrix(pbmc@raw.data[,which(colnames(pbmc@raw.data) %in% colnames(pbmc@scale.data))])
pbmc.data=as.matrix(pbmc@scale.data)
used_gene=rownames(pbmc.data)

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
ONE=.data2one(pbmc.raw.data, used_gene, CPU=4, PCNUM=50, SEED=123,  PP=30)
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
saveRDS(MEAN,file='MEAN.RDS')
    
PMAT=getPMAT(EXP, LR, BIN, MEAN)
saveRDS(PMAT,file='PMAT.RDS')

CMAT=getCMAT(EXP,LR,PMAT,PRO=TRUE)
saveRDS(CMAT,file='CMAT.RDS')

pdf('2HEAT.pdf',width=15,height=13)
library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',
  col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
dev.off()

OUT=getPAIR(CMAT)
PAIR=OUT$PAIR
SCORE=OUT$SCORE
RANK=OUT$RANK
saveRDS(PAIR,file='PAIR.RDS')

#---- !!! Changed in Seurat 3.0 !!! ----
VEC=pbmc@dr$tsne@cell.embeddings
#--------------------------------------- 
# For Seurat 3.0, please use:
# VEC=pbmc@reductions$tsne@cell.embeddings
#---------------------------------------

pdf('3CPlot.pdf',width=12,height=10)
CPlot(VEC,PAIR[1:200,],BINTAG)
dev.off()


ORITAG=as.character(pbmc@ident)
ORITAG[which(ORITAG=='Pericyte/\nFibroblast')]='Pericyte Fibroblast'

NET=getNET(PAIR, BINTAG,ORITAG )
write.table(NET,file='NET.txt',sep='\t',row.names=F,col.names=T,quote=F)
   
CN=getCN(NET)
pdf('4DPlot.pdf',width=20,height=20)
DP=DPlot(NET, CN, COL=3)
dev.off()


SIG_INDEX=which(DP<0.05)
SIG_PAIR=names(SIG_INDEX)

pdf('5LPlot.pdf',width=20,height=20)
RCN=trunc(sqrt(length(SIG_PAIR))+1)
par(mfrow=c(RCN,RCN))
i=1
while(i<= length(SIG_PAIR) ){
    this_pair=SIG_PAIR[i]
    LT=unlist(strsplit(this_pair, "_to_"))[1]
    RT=unlist(strsplit(this_pair, "_to_"))[2]
    LP=LPlot(LT, RT, NET, PMAT,SEED=123)    
    colnames(LP)=paste0(c('Lexp','Rexp'),'_',c(LT,RT))
    write.table(LP,file=paste0(as.character(SIG_INDEX[i]),'.tsv'),row.names=T,col.names=T,sep='\t',quote=F)
    print(i)
    i=i+1}
dev.off()













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














