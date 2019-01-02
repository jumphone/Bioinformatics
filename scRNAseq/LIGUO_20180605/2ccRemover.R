#https://cran.r-project.org/web/packages/ccRemover/vignettes/ccRemover_tutorial.html

library(Seurat)
library(ccRemover)

exp_data=read.table('5000_run1663_normalized.txt.rmdup',header=T,row.names=1,sep='\t')

#### ccRemover #############
mean_gene_exp <- rowMeans(exp_data)
exp_data_cen <- exp_data - mean_gene_exp 
gene_names <- rownames(exp_data_cen)
cell_cycle_gene_indices <- gene_indexer(gene_names, species = "mouse", name_type = "symbols" )
length(cell_cycle_gene_indices)
#1377
if_cc <- rep(FALSE,nrow(exp_data_cen)) 
if_cc[cell_cycle_gene_indices] <- TRUE
summary(if_cc)
dat <- list(x=exp_data_cen, if_cc=if_cc)
xhat <- ccRemover(dat, bar=FALSE)
xhat <- xhat + mean_gene_exp
############################

#EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes=0)
EXP = CreateSeuratObject(raw.data = xhat, min.cells = 0, min.genes=0)



mito.genes <- grep(pattern = "^mt-", x = rownames(x = EXP@data), value = TRUE)
percent.mito <- colSums(EXP@data[mito.genes, ]) / colSums(EXP@data)
EXP <- AddMetaData(object = EXP, metadata = percent.mito, col.name = "percent.mito")


pdf('Seurat_QC.pdf',width=30,height=15)
VlnPlot(object = EXP, features.plot = c("nGene", "nUMI",'percent.mito'), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "nGene")
dev.off()

EXP=FilterCells(object = EXP, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(3000, 0.05))

length(EXP@data[1,])

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

pdf('Seurat_VarGene.pdf')
EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.5)
dev.off()

length(x=EXP@var.genes)

'Prom1' %in% rownames(EXP@data)
'Prom1' %in% EXP@var.genes #F
'Fut4' %in% EXP@var.genes
'Nes' %in% EXP@var.genes
'Egfr' %in% EXP@var.genes
'Cd15' %in% EXP@var.genes #F
'Slc1a3' %in% EXP@var.genes
'Sox2' %in% EXP@var.genes
'Fabp7' %in% EXP@var.genes
'Nr2e1' %in% EXP@var.genes
'Id3' %in% EXP@var.genes
'Clu' %in% EXP@var.genes
'Sox9' %in% EXP@var.genes #F
'Vcam1' %in% EXP@var.genes
'Slc1a2' %in% EXP@var.genes
'Id2' %in% EXP@var.genes
'Sox11' %in% EXP@var.genes
'Apoe' %in% EXP@var.genes
'Tbr2' %in% EXP@var.genes
'Ntsr2' %in% EXP@var.genes


#PrintPCA(object = EXP, pcs.print = 1:5, genes.print = 5, use.full = FALSE)


#Stem_gene=c('Prom1','Nes','Egfr','Cd15','Slc1a3','Sox2','Fabp7','Nr2e1','Id3','Clu','Sox9','Vcam1','Slc1a2','Id2','Sox11','Apoe','Tbr2','Ntsr2')

EXP = ScaleData(object = EXP,vars.to.regress = c("percent.mito", "nUMI"), genes.use = EXP@var.genes)

#EXP = ScaleData(object = EXP,genes.use = EXP@var.genes)

PCNUM=40
EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )


PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = EXP, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
#EXP = ScaleData(object = EXP,vars.to.regress = c("percent.mito", "nUMI"), genes.use = Stem_gene)
#EXP <- RunPCA(object = EXP, pc.genes = Stem_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )



PCElbowPlot(object = EXP,num.pc=PCNUM)


PrintPCA(object = EXP, pcs.print = 1:20)



PCUSE=1:6
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE, check_duplicates = FALSE)

RES=0.5
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)

FeaturePlot(object = EXP, features.plot = c('Sox2','Egfr','Id3','Olig2','Zic1','Apoe','Gfap','Slc1a3','Nes'), cols.use = c("grey", "red"), reduction.use = "tsne")
FeaturePlot(object = EXP, features.plot = c('Sox2'),cols.use = c("grey", "red"), reduction.use = "tsne")


#FeaturePlot(object = EXP, features.plot = c('Cd133'),cols.use = c("grey", "red"), reduction.use = "tsne")

#save(EXP, file = "EXP.Robj")
load('EXP.Robj')

pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)


pdf('Result.pdf',width=20,height=15)

TSNEPlot(object = EXP,do.label=T)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()




