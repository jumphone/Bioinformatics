
library(Seurat)
pbmc.data=read.delim('matt_combined.txt',sep='\t',row.names=1,header=T)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")

VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI"), nCol = 2)

pbmc <- FilterCells(object = pbmc, subset.names = c("nUMI"), low.thresholds = c(1000000))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)



pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc)

PCNUM=10
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, pcs.compute=PCNUM, do.print = TRUE, pcs.print = 1:5, genes.print = 5)


PCUSE=1:10
pbmc <- RunTSNE(object = pbmc, perplexity=5, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE)


pbmc@meta.data$tag=colnames(pbmc@data)
TSNEPlot(object = pbmc,pt.size=4, do.label=T)
TSNEPlot(object = pbmc,pt.size=4,group.by='tag', do.label=T)
#PCAPlot(object = pbmc,pt.size=4, do.label=T)
