library(Seurat)
library(dplyr)

pbmc.data <- Read10X(data.dir = "./hg19/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 100,  project = "DIPG")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")


pdf('QC.pdf',width=10,height=6)
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
dev.off()



