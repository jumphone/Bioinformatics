library(Seurat)
library(dplyr)

case.data <- read.table('CASE.txt',sep='\t',check.name=F,row.names=1,header=T)
case <- CreateSeuratObject(raw.data = case.data, min.cells = 3, min.genes = 200, project = "Natalie")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = case@data), value = TRUE)
percent.mito <- Matrix::colSums(case@raw.data[mito.genes, ])/Matrix::colSums(case@raw.data)
case <- AddMetaData(object = case, metadata = percent.mito, col.name = "percent.mito")
case <- AddMetaData(object = case, metadata = case@ident, col.name = "batch")
case@meta.data$stim <- "case"
case=FilterCells(object = case, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(4000, 0.2))
case <- NormalizeData(object = case, normalization.method = "LogNormalize", scale.factor = 10000)
case = ScaleData(object = case,vars.to.regress = c("percent.mito", "nUMI", "batch"))

wt.data <- read.table('WT.txt',sep='\t',check.name=F,row.names=1,header=T)
wt <- CreateSeuratObject(raw.data = wt.data, min.cells = 3, min.genes = 200, project = "Natalie")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = wt@data), value = TRUE)
percent.mito <- Matrix::colSums(wt@raw.data[mito.genes, ])/Matrix::colSums(wt@raw.data)
wt <- AddMetaData(object = wt, metadata = percent.mito, col.name = "percent.mito")
wt <- AddMetaData(object = wt, metadata = wt@ident, col.name = "batch")
wt@meta.data$stim <- "wt"
wt=FilterCells(object = wt, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(4000, 0.2))
wt <- NormalizeData(object = wt, normalization.method = "LogNormalize", scale.factor = 10000)
wt = ScaleData(object = wt,vars.to.regress = c("percent.mito", "nUMI", "batch"))

case <- FindVariableGenes(case, do.plot = F)
wt <- FindVariableGenes(wt, do.plot = F)

g.1 <- head(rownames(case@hvg.info), 1500)
g.2 <- head(rownames(wt@hvg.info), 1500)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(case@scale.data))
genes.use <- intersect(genes.use, rownames(wt@scale.data))



combined_data = RunCCA(case, wt, genes.use = genes.use, num.cc = 30)




PCNUM=40
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

PCUSE=1:35
pbmc=RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)

pdf('./TSNE.pdf',width=30,height=15)
TSNEPlot(object = pbmc)
dev.off()

save(pbmc, file = "./pbmc.Robj")


