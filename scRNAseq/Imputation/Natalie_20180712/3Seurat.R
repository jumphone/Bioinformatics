library(Seurat)
library(dplyr)
library(Matrix)

exp_data=read.table('magic.csv',header=T,row.names=1,sep='\t',check.names=F)
exp_data=t(exp_data)

pbmc <- CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes = 0, project = "Project")

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

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.3, print.output = 0, save.SNN = TRUE,force.recalc=T)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

pdf('TSNE.pdf')
TSNEPlot(object = pbmc,do.label = TRUE)
dev.off()


#save(pbmc,file='SAVE.Robj')

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,test.use='t', logfc.threshold=0.01)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf('HEAT.pdf',width=15,height=15)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
write.table(top10,file='top10.txt',sep='\t',quote=F,row.names=T,col.names=T)
write.table(pbmc.markers,file='markers.txt',sep='\t',quote=F,row.names=T,col.names=T)
#####################

ori_exp_data=read.table('MATRIX.txt',header=T,row.names=1,sep=',',check.names=F)
ori_exp_data=t(ori_exp_data)

ori_pbmc <- CreateSeuratObject(raw.data = ori_exp_data, min.cells = 0, min.genes = 0, project = "Origin")
ori_pbmc <- NormalizeData(object = ori_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
ori_pbmc <- ScaleData(object = ori_pbmc)
ori_pbmc@ident=pbmc@ident

pbmc.markers <- FindAllMarkers(object = ori_pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,test.use='t', logfc.threshold=0.01)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf('ori_HEAT.pdf',width=15,height=15)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
write.table(top10,file='ori_top10.txt',sep='\t',quote=F,row.names=T,col.names=T)
write.table(pbmc.markers,file='ori_markers.txt',sep='\t',quote=F,row.names=T,col.names=T)


###########################




