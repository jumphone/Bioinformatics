library(Seurat)
library(cowplot)


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

D1 <- Read10X(data.dir = "./CDC42_HET/")
D2 <- Read10X(data.dir = "./Small_Intestine/")
D1=as.matrix(D1)
D2=as.matrix(D2)
colnames(D1)=paste0('HET_', colnames(D1))
colnames(D2)=paste0('KO_', colnames(D2))

DATA=.simple_combine(D1,D2)$combine
BATCH=rep('KO',ncol(DATA))
BATCH[c(1:ncol(D1))]='HET'




mybeer=BEER(DATA,BATCH,MTTAG='^mt-')

PCUSE=mybeer$select#.getUSE(mybeer, 0.75,0.75)#mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

pbmc=mybeer$seurat
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)

DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)    











ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)

DimPlot(immune.combined, reduction = "umap", group.by = "stim")

