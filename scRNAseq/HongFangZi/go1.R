

setwd('/users/zha8dh/tianlab/HFZ')
source('./BEER.R')
DATA=readRDS('./DATA.RDS')
BATCH=readRDS('./BATCH.RDS')

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   
saveRDS(mybeer,file='mybeer.RDS')


##################################



setwd('/users/zha8dh/tianlab/HFZ')
source('./BEER.R')
.set_python('~/anaconda3/bin/python3')
mybeer=readRDS('./mybeer.RDS')

pdf('~/Downloads/HFZ1.pdf',width=10,height=10)

PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

dev.off()

#############
pbmc <- mybeer$seurat
################

pdf('~/Downloads/HFZ2.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.5,label=F)
dev.off()


pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=20, NT=10)
pbmc@reductions$umap@cell.embeddings=umap

pdf('~/Downloads/HFZ3.10.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.5,label=F)
dev.off()


saveRDS(pbmc,file='pbmc.RDS')








