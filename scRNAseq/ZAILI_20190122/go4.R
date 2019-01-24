library(Seurat)
library(dplyr)

pbmc.data <- readRDS('cerebullum_dev.RDS')
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
      
length(x = pbmc@var.genes)
#4602

pbmc <- ScaleData(object = pbmc, genes.use =pbmc@var.genes, vars.to.regress = c("nUMI", "percent.mito"),num.cores =4,do.par=TRUE)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

meta.data=read.table('cerebellum_cell_metadata.tsv',sep='\t',header=T)

