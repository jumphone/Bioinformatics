system('curl https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R > scRef.R')

library(Seurat)
source('scRef')

tag_data=read.table('GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt',sep='\t',row.names=1,header=T)


pbmc.data=read.table('GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt',sep='\t',row.names=1,header=T)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

PCNUM=50
pbmc <- RunPCA(object = pbmc,pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)

pbmc@meta.data$paper.type=tag_data[,6]
pbmc@meta.data$paper.cluster=tag_data[,5]
pbmc@meta.data$paper.condition=tag_data[,3]


PCUSE=1:20
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE,do.fast=T)



DimPlot(pbmc, group.by='orig.ident', reduction.use='tsne')
DimPlot(pbmc, group.by='paper.type', reduction.use='tsne')
DimPlot(pbmc, group.by='paper.cluster', reduction.use='tsne')
DimPlot(pbmc, group.by='paper.condition', reduction.use='tsne')





#pbmc <- RunUMAP(object = pbmc, dims.use = PCUSE)
#DimPlot(pbmc, group.by='orig.ident', reduction.use='umap')
#pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = PCUSE, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)



