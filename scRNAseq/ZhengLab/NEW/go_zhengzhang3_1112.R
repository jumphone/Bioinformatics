

setwd('F:/Zhenglab/NewZhengZhang3')
library(Seurat)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#REF=readRDS('REF.RDS')
#LABEL=readRDS('LABEL.RDS')



CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO <- Read10X(data.dir = "./CDC42_KO")
CDC42Rescue <- Read10X(data.dir = "./CDC42Rescue")


BATCH=c( rep('CDC42HET',ncol(CDC42HET )), rep('CDC42KO',ncol(CDC42KO )),
        rep('CDC42Rescue',ncol(CDC42Rescue )))


D1=.simple_combine(CDC42HET, CDC42KO)$combine

DATA=.simple_combine(D1, CDC42Rescue)$combine


#############################

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE )

# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

pbmc <- mybeer$seurat
PCUSE=mybeer$select   
#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)

saveRDS(mybeer,file='mybeer.RDS')
saveRDS(pbmc,file='pbmc.RDS')


##############
DimPlot(pbmc, reduction = "umap", split.by = "batch",ncol=2)

########################

VEC=pbmc@reductions$umap@cell.embeddings

# Here, we use K-means to do the clustering
N=200
set.seed(123)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)+ NoLegend()

pbmc@meta.data$celltype=rep('Enterocyte',ncol(pbmc))

##################################################









