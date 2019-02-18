
library(Seurat)
library(dplyr)


pbmc.data=read.table('decidua-20181218_exon_tagged.dge.txt',sep='\t',header=T,row.names=1)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.2))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

PCNUM=50
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, pcs.compute=PCNUM, genes.print = 5)

PCElbowPlot(object = pbmc,num.pc=PCNUM)

PCUSE=1:6
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pbmc <- RunUMAP(pbmc, dims.use = PCUSE)

TSNEPlot(object = pbmc)
#DimPlot(pbmc, reduction.use = "umap")

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = PCUSE, resolution = 0.6, print.output = 0, save.SNN = TRUE)

TSNEPlot(object = pbmc)
#DimPlot(pbmc, reduction.use = "umap")

saveRDS(pbmc,file='orig_decidua.RDS')

#########################
source


