library(Seurat)
library(dplyr)

pbmc.data <- read.table('ALL.txt',sep='\t',check.name=F,row.names=1,header=T)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "Natalie")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pbmc <- AddMetaData(object = pbmc, metadata = pbmc@ident, col.name = "batch")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

pbmc=FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(2500, 0.1))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot = F, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.5)
length(x=pbmc@var.genes) #10963
pbmc = ScaleData(object = pbmc,vars.to.regress = c("percent.mito", "nUMI", "batch"), genes.use=pbmc@var.genes)

PCNUM=40
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.compute=PCNUM, pcs.print = 1:5,  genes.print = 5)




