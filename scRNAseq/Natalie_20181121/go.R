library('Seurat')

##########################################
GBM.data=read.table('./GBM/GBM_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
GBM <- CreateSeuratObject(raw.data = GBM.data, min.cells = 3, min.genes = 10, project = "GBM")

##########################################
medullo.data=read.table('./medullo/medullo_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
medullo <- CreateSeuratObject(raw.data = medullo.data, min.cells = 3, min.genes = 10, project = "medullo")

##########################################
nb6.data=read.table('./nb6/nb6_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
nb6 <- CreateSeuratObject(raw.data = nb6.data, min.cells = 3, min.genes = 10, project = "nb6")

##########################################
rm(GBM.data)
rm(medullo.data)
rm(nb6.data)
save.image('raw.RData')
##########################################

pdf('QC.pdf')
VlnPlot(object = GBM, features.plot = c("nGene", "nUMI"), nCol = 2)
VlnPlot(object = medullo, features.plot = c("nGene", "nUMI"), nCol = 2)
VlnPlot(object = nb6, features.plot = c("nGene", "nUMI"), nCol = 2)
dev.off()

GBM <- FilterCells(object = GBM, subset.names = c("nGene"),  low.thresholds = c(200), high.thresholds = c(2500))
medullo <- FilterCells(object = medullo, subset.names = c("nGene"),  low.thresholds = c(50), high.thresholds = c(1000))
nb6 <- FilterCells(object = nb6, subset.names = c("nGene"),  low.thresholds = c(50), high.thresholds = c(1000))

GBM <- NormalizeData(object = GBM, normalization.method = "LogNormalize", scale.factor = 10000)
medullo <- NormalizeData(object = medullo, normalization.method = "LogNormalize", scale.factor = 10000)
nb6 <- NormalizeData(object = nb6, normalization.method = "LogNormalize", scale.factor = 10000)

GBM <- FindVariableGenes(object = GBM, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
medullo <- FindVariableGenes(object = medullo, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
nb6  <- FindVariableGenes(object = nb6 , mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)

length(x = GBM@var.genes)
length(x = medullo@var.genes)
length(x = nb6@var.genes)

GBM <- ScaleData(object = GBM, vars.to.regress = c("nUMI"))
medullo <- ScaleData(object = medullo, vars.to.regress = c("nUMI"))
nb6 <- ScaleData(object = nb6, vars.to.regress = c("nUMI"))














