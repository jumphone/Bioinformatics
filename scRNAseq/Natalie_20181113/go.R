library(Seurat)
library(dplyr)

pbmc.data <- read.table('EXP12345678.combined.txt',sep='\t',check.name=F,row.names=1,header=T)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "Natalie")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pbmc <- AddMetaData(object = pbmc, metadata = pbmc@ident, col.name = "batch")

pbmc=FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(4000, 0.1))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pdf('./VarGene.pdf')
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.8)
dev.off()

length(x=pbmc@var.genes)
#2975

pbmc = ScaleData(object = pbmc,vars.to.regress = c("percent.mito", "nUMI", "batch"), genes.use = pbmc@var.genes)



