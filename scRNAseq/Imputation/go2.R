
library(Seurat)
library(dplyr)
library(Matrix)

exp_data=read.table('magic.csv',header=T,row.names=1,sep='\t',check.names=F)
exp_data=t(exp_data)

pbmc <- CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes = 0, project = "DIPG")


mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

pdf('QC.pdf',width=10,height=6)
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
dev.off()

pdf('VAR.pdf')
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

pdf('PCA.pdf')
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
PCAPlot(object = pbmc, dim.1 = 2, dim.2 = 3)
PCAPlot(object = pbmc, dim.1 = 3, dim.2 = 4)
PCElbowPlot(object = pbmc)
dev.off()

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:4, resolution = 0.2, print.output = 0, save.SNN = TRUE)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:4, do.fast = TRUE)

pdf('TSNE.pdf')
TSNEPlot(object = pbmc,do.label = TRUE)
dev.off()

