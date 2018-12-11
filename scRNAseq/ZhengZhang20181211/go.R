library(Seurat)
library(dplyr)

source('scRef.R')
pbmc.data <- Read10X(data.dir = "./WT/mm10/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pdf('WT_QC.pdf')
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(4000, 0.5))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 2, y.cutoff = 0.25)
length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

PCNUM=150
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes,pcs.compute = PCNUM, do.print = F, pcs.print = 1:5, 
    genes.print = 5)

pdf('WT_PCA.pdf')
PCElbowPlot(object = pbmc)
dev.off()

PCUSE=1:25
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
TSNEPlot(object = pbmc)

pdf('WT_TSNE.pdf')
TSNEPlot(object = pbmc)
dev.off()

saveRDS(pbmc,file='WT.RDS')

