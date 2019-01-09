library(Seurat)
library(dplyr)
pbmc.data <- Read10X(data.dir='./10x\ data/hg19/')
pbmc = CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200,project = "10X_PBMC")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
mito.genes
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05)) 
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
pbmc <- ScaleData(object = pbmc,genes.use = pbmc@var.genes)
length(x = pbmc@var.genes)

EXP=pbmc

################
library(Seurat)
library(ccRemover)
var_exp_data=as.matrix(EXP@data[which(rownames(EXP@data) %in% EXP@var.genes),])
mean_gene_exp <- rowMeans(var_exp_data)
exp_data_cen <- var_exp_data - mean_gene_exp 
gene_names <- rownames(exp_data_cen)
cell_cycle_gene_indices <- gene_indexer(gene_names, species = "mouse", name_type = "symbols" )
length(cell_cycle_gene_indices)
#1377
if_cc <- rep(FALSE,nrow(exp_data_cen)) 
if_cc[cell_cycle_gene_indices] <- TRUE
summary(if_cc)
exp_data_cen=as.matrix(exp_data_cen)
dat <- list(x=exp_data_cen, if_cc=if_cc)
xhat <- ccRemover(dat, bar=FALSE)
xhat <- xhat + mean_gene_exp
saveRDS(xhat,file='xhat.RDS')
#########

xhat=readRDS('xhat.RDS')
deNeg<-function(X){
   return(X-min(X))
}
pxhat=t(apply(xhat,1,deNeg))
colnames(pxhat)=colnames(xhat)
rownames(pxhat)=rownames(xhat)
pbmc = CreateSeuratObject(raw.data = pxhat, min.cells = 0, min.genes=0)
pbmc = NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc = ScaleData(object = pbmc)
EXP@scale.data=pbmc@scale.data

#########

pbmc = EXP

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, pcs.compute=100, genes.print = 5)
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
PCElbowPlot(object = pbmc)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:80, resolution = 0.6, print.output = 0, save.SNN = TRUE) 
PrintFindClustersParams(object = pbmc) 
pbmc <- RunTSNE(object = pbmc, dims.use = 1:80, do.fast = TRUE)
TSNEPlot(object = pbmc)
save.image('DATA.RData')

