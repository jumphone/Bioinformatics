library(Seurat)
library(dplyr)

# Load the PBMC dataset
require(data.table)
pbmc.data=data.frame(fread("E6701_exp.txt",header=T), row.names=1)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 100, min.genes = 1000, project = "10X_PBMC")

rm(pbmc.data)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(5000, 0.1))

#V=apply(pbmc@data,1,var)

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F, mean.function = ExpMean, dispersion.function = LogVMR,  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI", "percent.mito"))



