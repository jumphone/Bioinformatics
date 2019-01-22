library(Seurat)
source('scRef.R')

load('R4.RObj')
PCUSE=1:15
pbmc=RunUMAP(pbmc, reduction.use = "pca", dims.use = PCUSE)
DimPlot(pbmc, reduction.use = "umap")
saveRDS(pbmc,file='R4_umap.RDS')

load('R18058913.RObj')
PCUSE=1:15
pbmc=RunUMAP(pbmc, reduction.use = "pca", dims.use = PCUSE)
DimPlot(pbmc, reduction.use = "umap")
saveRDS(pbmc,file='R18058913_umap.RDS')

load('R18058914.RObj')
PCUSE=1:15
pbmc=RunUMAP(pbmc, reduction.use = "pca", dims.use = PCUSE)
DimPlot(pbmc, reduction.use = "umap")
saveRDS(pbmc,file='R18058914_umap.RDS')

load('R18059655.RObj')
PCUSE=1:15
pbmc=RunUMAP(pbmc, reduction.use = "pca", dims.use = PCUSE)
DimPlot(pbmc, reduction.use = "umap")
saveRDS(pbmc,file='R18059655_umap.RDS')
