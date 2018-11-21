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

PCNUM=20
GBM <- RunPCA(object = GBM, pcs.compute=PCNUM, pc.genes = GBM@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
medullo <- RunPCA(object = medullo, pcs.compute=PCNUM, pc.genes = medullo@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
nb6 <- RunPCA(object = nb6, pcs.compute=PCNUM, pc.genes = nb6@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)

pdf('PCA.pdf')
PCElbowPlot(object = GBM)
PCElbowPlot(object = medullo)
PCElbowPlot(object = nb6)
dev.off()

PCUSE=1:20
GBM <- FindClusters(object = GBM, reduction.type = "pca", dims.use = PCUSE, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)
medullo <- FindClusters(object = medullo, reduction.type = "pca", dims.use = PCUSE, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)
nb6 <- FindClusters(object = nb6, reduction.type = "pca", dims.use = PCUSE, 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)

GBM <- RunTSNE(object = GBM, dims.use = PCUSE, do.fast = TRUE)
medullo  <- RunTSNE(object = medullo, dims.use = PCUSE, do.fast = TRUE)
nb6 <- RunTSNE(object = nb6, dims.use = PCUSE, do.fast = TRUE)

pdf('TSNE.pdf')
TSNEPlot(object = GBM)
TSNEPlot(object = medullo)
TSNEPlot(object = nb6)
dev.off()

GBM.markers <- FindAllMarkers(object = GBM, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
medullo.markers <- FindAllMarkers(object = medullo, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
nb6.markers <- FindAllMarkers(object = nb6, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

library(dplyr)
GBM.top10 <- GBM.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
medullo.top10 <- medullo.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
nb6.top10 <- nb6.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

pdf('HEAT.pdf',width=10,height=12)
DoHeatmap(object = GBM, genes.use = GBM.top10$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = medullo, genes.use = medullo.top10$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = nb6, genes.use = nb6.top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

save.image('done.RData')
