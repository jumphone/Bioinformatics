

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
dim(mybeer$seurat)
#13828 64000
################################################

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")



pdf('~/Downloads/HFZ_QC.pdf',width=12,height=7)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

dim(pbmc)
#[1] 13828 36797
saveRDS(as.matrix(pbmc@assays$RNA@counts),file='DATA.QC.RDS')
saveRDS(pbmc@meta.data$batch,file='BATCH.QC.RDS')

################################################################

DATA.QC=as.matrix(pbmc@assays$RNA@counts)
BATCH.QC=as.matrix(pbmc@meta.data$batch)

rm(mybeer)
rm(pbmc)
gc()

mybeer.qc=BEER(DATA.QC, BATCH.QC, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   
saveRDS(mybeer.qc,file='mybeer.qc.RDS')




pbmc <- mybeer.qc$seurat
PCUSE=mybeer.qc$select

pdf('~/Downloads/HFZ2.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.5,label=F)
dev.off()


pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
pbmc <- RunUMAP(pbmc, dims = PCUSE)

saveRDS(pbmc,file='pbmc.RDS')


##umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
#pbmc@reductions$umap@cell.embeddings=umap

pdf('~/Downloads/HFZ3.combat.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.5,label=F)
dev.off()


pdf('~/Downloads/HFZ4.EXP.pdf',width=10,height=10)
FeaturePlot(pbmc, features=c('PTPRC','YAP1') )
dev.off()


VEC=pbmc@reductions$umap@cell.embeddings

# Here, we use K-means to do the clustering
N=50
set.seed(1)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)

pdf('~/Downloads/HFZ5.CLUST.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)
dev.off()



Idents(pbmc)=pbmc@meta.data$clust
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


print('ok')






