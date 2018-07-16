library(Seurat)
library(dplyr)

GENE_CUTOFF_LOW = 100
GENE_CUTOFF_HIGH = 1000
CELL_CUTOFF_LOW = 10
MITO_CUTOFF_HIGH = 0.1


pbmc.data <- Read10X(data.dir = "./raw_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = CELL_CUTOFF_LOW, min.genes = GENE_CUTOFF_LOW,  project = "Human")
dim(pbmc@raw.data)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")


pdf('QC.pdf',width=10,height=6)
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
dev.off()

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(GENE_CUTOFF_LOW, -Inf), high.thresholds = c(GENE_CUTOFF_HIGH, MITO_CUTOFF_HIGH))

dim(pbmc@data)


save(pbmc, file = "Seurat.Robj")
output = t(as.matrix(pbmc@data))
write.table(output,file='MATRIX.txt',quote=F,sep=',',row.names=T,col.names=T)

