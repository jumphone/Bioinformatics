
setwd('/Volumes/Feng/Zhenglab/Combine')


setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
.set_python('C:/Users/cchmc/Anaconda3/python.exe')

DATA=readRDS('DATA2.RDS') # Kodanda
BATCH=readRDS('BATCH2.RDS')



mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=3, GN=5000, SEED=1, COMBAT=TRUE)
saveRDS(mybeer,file='mybeer2.RDS')




PCUSE=mybeer$select

COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


pbmc <- mybeer$seurat  
PCUSE=mybeer$select
#PCUSE=.selectUSE(mybeer, CUTR=0.8, CUTL=0.8, RR=0.5, RL=0.5)

#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)

saveRDS(pbmc,file='pbmc2.RDS')

















