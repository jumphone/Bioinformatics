library(Seurat)

exp_data=read.table('ctrl.csv',sep=',',header=T,row.names=1)

dim(exp_data)

EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=200)

dim(EXP@raw.data)




pdf('Seurat_QC.pdf',width=30,height=15)
VlnPlot(object = EXP, features.plot = c("nGene", "nUMI"), nCol = 2)
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "nGene")
dev.off()


EXP=FilterCells(object = EXP, subset.names = c("nGene"), low.thresholds = c(500), high.thresholds = c(4000))


dim(EXP@data)

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

pdf('Seurat_VarGene.pdf')
EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.5)
dev.off()

length(x=EXP@var.genes)

EXP = ScaleData(object = EXP,vars.to.regress = c( "nUMI"), genes.use = EXP@var.genes)


PCNUM=40
EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )


PCUSE=1:6
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE, check_duplicates = FALSE)


RES=0.1
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)



pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

write.table(top10,'top10.txt',sep='\t',quote=F,row.names=T,col.names=T)
write.table(pbmc.markers,'markers.txt',sep='\t',quote=F,row.names=T,col.names=T)

pdf('Result.pdf',width=20,height=15)

TSNEPlot(object = EXP,do.label=T)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()



