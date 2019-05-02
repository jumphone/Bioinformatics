library('Seurat')
library(dplyr)
library(Matrix)

source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')


load('Seurat_EXP_cluster.Robj')

pbmc=EXP_cluster
used=which(as.numeric(as.character(EXP_cluster@ident)) %in% c(2,9,14,17,19,23))

pbmc.raw.data=getSeuratRAW(pbmc@raw.data,pbmc@scale.data)[,used]

############
tmp=CreateSeuratObject(raw.data = pbmc.raw.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")
mito.genes <- grep(pattern = "^Mt", x = rownames(x = tmp@data), value = TRUE)
percent.mito <- colSums(tmp@data[mito.genes, ]) / colSums(tmp@data)
tmp <- AddMetaData(object = tmp, metadata = percent.mito, col.name = "percent.mito")
tmp <- NormalizeData(object = tmp, normalization.method = "LogNormalize",  scale.factor = 10000)
tmp <- FindVariableGenes(object = tmp,do.plot=FALSE, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = tmp@var.genes)
tmp <- ScaleData(object = tmp, vars.to.regress = c("nUMI",'percent.mito'), genes.use = tmp@var.genes)
#################

pbmc.data=as.matrix(tmp@scale.data)

used_gene=tmp@var.genes
pbmc.raw.data=pbmc.raw.data[which(rownames(pbmc.raw.data) %in% used_gene),]
pbmc.data=pbmc.data[which(rownames(pbmc.data) %in% used_gene),]

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
ONE=.data2one(pbmc.raw.data, used_gene, CPU=4, PCNUM=50, SEED=123,  PP=30)
saveRDS(ONE,file='ONE.RDS')

OUT=getBIN(ONE)
BIN=OUT$BIN
BINTAG=OUT$TAG
saveRDS(BIN,file='BIN.RDS')
saveRDS(BINTAG,file='BINTAG.RDS')

pbmc@meta.data$bin=BINTAG
pdf('1ID.pdf',width=12,height=10)
DimPlot(pbmc,group.by='bin',reduction.use='tsne',do.label=T)
dev.off()

LR=read.table('RL_mouse.txt',header=T,sep='\t')
EXP=pbmc.data

MEAN=getMEAN(EXP, LR)
saveRDS(MEAN,file='MEAN.RDS')
    
PMAT=getPMAT(EXP, LR, BIN, MEAN)
saveRDS(PMAT,file='PMAT.RDS')

CMAT=getCMAT(EXP,LR,PMAT,PLUS=TRUE)
saveRDS(CMAT,file='CMAT.RDS')

pdf('2HEAT.pdf',width=12,height=10)
library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',
  col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
dev.off()

OUT=getPAIR(CMAT)
PAIR=OUT$PAIR
SCORE=OUT$SCORE
RANK=OUT$RANK
saveRDS(PAIR,file='PAIR.RDS')

VEC=pbmc@dr$tsne@cell.embeddings
pdf('3CPlot.pdf',width=12,height=10)
CPlot(VEC,PAIR[1:200,],BINTAG)
dev.off()

ORITAG=as.character(pbmc@ident)
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
    LP=LPlot(LT, RT, NET, PMAT,MAIN=as.character(SIG_INDEX[i]),SEED=123)    
    colnames(LP)=paste0(c('Lexp','Rexp'),'_',c(LT,RT))
    write.table(LP,file=paste0(as.character(SIG_INDEX[i]),'.tsv'),row.names=T,col.names=T,sep='\t',quote=F)
    print(i)
    i=i+1}
dev.off()




