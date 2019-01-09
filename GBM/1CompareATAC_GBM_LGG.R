TAG=as.character(read.table('GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.TAG.txt',header=F)[,1])
#a=read.table('GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.txt',row.names=1,header=T,sep='\t',check.names = F)
a=readRDS('GBMLGG_deNeg_ATAC_PanCan_Log2Norm_Counts.RDS')

library(Seurat)
library(dplyr)

pbmc <- CreateSeuratObject(raw.data = a, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc@meta.data$tag=TAG
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, do.plot=F,
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc)


PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)

PCElbowPlot(object = pbmc)

PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, perplexity=5, dims.use = PCUSE, do.fast = TRUE)

TSNEPlot(object = pbmc,group.by='tag')
