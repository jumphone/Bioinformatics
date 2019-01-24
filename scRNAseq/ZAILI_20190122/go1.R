
library('Seurat')
source('scRef.R')


R4.data <- Read10X(data.dir = "Group_R4/filtered_gene_bc_matrices/GRCh38")


pbmc <- CreateSeuratObject(raw.data = R4.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute = PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc)
PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pdf('R4_TSNE.pdf')
TSNEPlot(pbmc,no.legend=T)
dev.off()
save(file='R4.RObj',pbmc)

################################################################################



R18058913.data <- Read10X(data.dir = "Group_R18058913/filtered_gene_bc_matrices/hg19")


pbmc <- CreateSeuratObject(raw.data = R18058913.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute = PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc)
PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pdf('R18058913_TSNE.pdf')
TSNEPlot(pbmc,no.legend=T)
dev.off()
save(file='R18058913.RObj',pbmc)



################################################################################

R18059655.data <- Read10X(data.dir = "Group_R18059655/filtered_gene_bc_matrices/hg19")


pbmc <- CreateSeuratObject(raw.data = R18059655.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute = PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc)
PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pdf('R18059655_TSNE.pdf')
TSNEPlot(pbmc,no.legend=T)
dev.off()
save(file='R18059655.RObj',pbmc)

################################################################################

R18058914.data <- Read10X(data.dir = "Group_R18058914/filtered_gene_bc_matrices/hg19")

pbmc <- CreateSeuratObject(raw.data = R18058914.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute = PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc)
PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pdf('R18058914_TSNE.pdf')
TSNEPlot(pbmc,no.legend=T)
dev.off()
save(file='R18058914.RObj',pbmc)



################################################################################
################################################################################

BT311.data <- Read10X(data.dir = "BT311_filtered_feature_bc_matrix/filtered_feature_bc_matrix/GRCh38")

pbmc <- CreateSeuratObject(raw.data = BT311.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute = PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc)
PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pbmc=RunUMAP(pbmc, reduction.use = "pca", dims.use = PCUSE)

saveRDS(pbmc,file='BT311_umap.RDS')



################################################################################
################################################################################

BT310.data <- Read10X(data.dir = "BT310 _filtered_gene_bc_matrices/filtered_feature_bc_matrix/GRCh38")

pbmc <- CreateSeuratObject(raw.data = BT310.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute = PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc)
PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pbmc=RunUMAP(pbmc, reduction.use = "pca", dims.use = PCUSE)

saveRDS(pbmc,file='BT310_umap.RDS')

