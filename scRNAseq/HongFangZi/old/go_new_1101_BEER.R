setwd('/data/tianlab/zhangfeng/HFZ')
print(getwd())

source('BEER.R')
#####################
DATA=readRDS('./DATA.RDS')
BATCH=readRDS('./BATCH.RDS')


pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
Idents(pbmc)=BATCH
pbmc@meta.data$batch=BATCH
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


pdf('~/Downloads/HFZ_QC.pdf',width=12,height=7)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)
#########


BATCH=pbmc@meta.data$batch
DATA=as.matrix(pbmc@assays$RNA@counts[,which(colnames(pbmc@assays$RNA@counts) %in% colnames(pbmc@assays$RNA@data))])

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)
saveRDS(mybeer,file='mybeer.RDS')




pdf('~/Downloads/HFZ1.pdf',width=10,height=10)

PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

dev.off()

pbmc <- mybeer$seurat
PCUSE=mybeer$select
#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=10, NT=10)
pbmc@reductions$umap@cell.embeddings=umap

pdf('~/Downloads/HFZ2.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()

saveRDS(pbmc, file='pbmc.enh.RDS')



