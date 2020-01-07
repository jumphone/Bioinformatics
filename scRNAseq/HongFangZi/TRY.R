
setwd('F:/Hongfangzi')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

DATA=readRDS( file='./DATA.RDS')
BATCH=readRDS( file='./BATCH.RDS')

USED=which(BATCH %in% c('placenta0423.2','placenta508.1',
                       'placenta508.2','placenta514.2',
                       'placenta2018c', 'placenta2019c',
                       'placenta20190402'
                       ))




DATA=DATA[,USED]
BATCH=BATCH[USED]

library(Seurat)
# DATA: Expression matrix. Rownames are gene names. Colnames are cell names.
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,'pbmc_placenta.RDS')


VEC = pbmc@reductions$umap@cell.embeddings
rownames(VEC) = colnames(pbmc)
PCA = pbmc@reductions$pca@cell.embeddings

Sys.time()
PCA=vector.rankPCA(PCA)
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
Sys.time()




FeaturePlot(pbmc,features='ACTA2')
################################################################################



setwd('F:/Hongfangzi')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

DATA=readRDS( file='./DATA.RDS')
BATCH=readRDS( file='./BATCH.RDS')

USED=which(BATCH %in% c('decidua0117','decidua0417.2',
                       'decidua2018c','decidua20190215',
                       'decidua20190420', 'decidua2019c',
                       'decidua508','decidua510'
                       ))




DATA=DATA[,USED]
BATCH=BATCH[USED]

library(Seurat)
# DATA: Expression matrix. Rownames are gene names. Colnames are cell names.
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,'pbmc_placenta.RDS')


VEC = pbmc@reductions$umap@cell.embeddings
rownames(VEC) = colnames(pbmc)
PCA = pbmc@reductions$pca@cell.embeddings






# FB : ACTA2, COL1A1






