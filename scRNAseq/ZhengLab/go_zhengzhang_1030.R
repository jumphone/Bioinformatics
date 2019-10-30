
setwd('F:/Zhenglab/NewZhengZhang')
library(Seurat)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

REF=readRDS('REF.RDS')
LABEL=readRDS('LABEL.RDS')



CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO <- Read10X(data.dir = "./CDC42_KO")
CDC42Rescue <- Read10X(data.dir = "./CDC42Rescue")
YapHet <- Read10X(data.dir = "./YapHet")


BATCH=c(rep('Nature',ncol(REF)), rep('CDC42HET',ncol(CDC42HET )), rep('CDC42KO',ncol(CDC42KO )),
        rep('CDC42Rescue',ncol(CDC42Rescue )), rep('YapHet',ncol(YapHet )))



D1=.simple_combine(REF, CDC42HET)$combine
D2=.simple_combine(CDC42KO, CDC42Rescue)$combine
D3=.simple_combine(D1, D2)$combine
DATA=.simple_combine(D3, YapHet)$combine


TAG=c(LABEL, rep(NA, (ncol(DATA)-ncol(REF))))

saveRDS(DATA,file='DATA.RDS')
saveRDS(BATCH,file='BATCH.RDS')
saveRDS(TAG,file='TAG.RDS')
#################################################################################################




setwd('F:/Zhenglab/NewZhengZhang')
library(Seurat)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')



DATA=readRDS(file='DATA.RDS')
BATCH=readRDS(file='BATCH.RDS')
TAG=readRDS(file='TAG.RDS')



pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
Idents(pbmc)=BATCH
pbmc@meta.data$batch=BATCH
pbmc@meta.data$tag=TAG
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)


BATCH=pbmc@meta.data$batch
TAG=pbmc@meta.data$tag
DATA=as.matrix(pbmc@assays$RNA@counts[,which(colnames(pbmc@assays$RNA@counts) %in% colnames(pbmc@assays$RNA@data))])



mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE )
mybeer$seurat@meta.data$tag=TAG
saveRDS(mybeer,file='mybeer.RDS')


mybeer <- ReBEER(mybeer, GNUM=30, PCNUM=150, ROUND=1, SEED=1, RMG=NULL) 

saveRDS(mybeer,file='mybeer150.RDS')


######################

# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
#######################


pbmc <- mybeer$seurat
PCUSE <- mybeer$select
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)

DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)  



pbmc <- mybeer$seurat
PCUSE=mybeer$select   
#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)


DimPlot(pbmc, reduction.use='umap', group.by='tag', pt.size=0.1,label=T)


saveRDS(pbmc,file='pbmc.RDS')






