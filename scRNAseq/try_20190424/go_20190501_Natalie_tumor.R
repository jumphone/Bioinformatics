source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')


library(Seurat)
library(dplyr)
library(Matrix)


load('Seurat_EXP_cluster.Robj')

#pbmc.raw.data=as.matrix(EXP_cluster@raw.data[,which(colnames(EXP_cluster@raw.data) %in% colnames(EXP_cluster@scale.data))])
#pbmc.data=as.matrix(EXP_cluster@scale.data)

#source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
#ONE=.data2one(pbmc.raw.data, rownames(pbmc.data), CPU=4, PCNUM=50, SEED=123,  PP=30)
ONE=readRDS('ONE.RDS')


BIN=readRDS('BIN.RDS')


pbmc=EXP_cluster
used=which(as.numeric(as.character(EXP_cluster@ident)) %in% c(2,9,14,17,19,23))



TAG=BIN2BINTAG(BIN,ONE)
BINTAG=rep(NA,length(pbmc@ident))
BINTAG[used]=TAG
#rep(NA,length(ONE))
#i=1
#while(i<=ncol(BIN)){
#BINTAG[BIN[,i]]=i
#i=i+1
#}


pbmc@meta.data$bin=BINTAG
png('ID.png',width=1200,height=1000)
DimPlot(pbmc,group.by='bin',reduction.use='tsne',do.label=T)
dev.off()



#EXP=pbmc.data
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
CPlot(VEC,PAIR[1:100,],BINTAG)







###################
LT='1.5mo_53,41,40_tumorSCs'
RT='30,32_Tcell'
ORITAG=rep('NA',length(pbmc@ident))
ORITAG[which(BINTAG %in% c(53,41,40))]=LT
ORITAG[which(BINTAG %in% c(30,32))]=RT

NET=getNET(PAIR[1:100,], BINTAG,ORITAG )
CN=getCN(NET)
DP=DPlot(NET, CN,COL=3,PLOT=F)
pdf(paste0(LT,'_to_',RT,'.pdf'),width=12,height=12)
LP=LPlot(LT,RT,NET,PMAT,SEED=12)
colnames(LP)=paste0(colnames(LP),'_',c(LT,RT))
write.table(LP,file=paste0(LT,'_',RT,'.tsv'),col.names=T,row.names=T,quote=F,sep='\t')
dev.off()
###################






