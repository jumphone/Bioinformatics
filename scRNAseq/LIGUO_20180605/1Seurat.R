library(Seurat)

PCNUM=20
PCUSE=1:15
RES=0.8

exp_data=read.table('run1642_10000.dge.txt',header=T,row.names=1)


#exp_data=read.table('N709_6000_picard.bam.clean.bam.dge.txt',header=T,row.names=1)


EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=200)

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
EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0.0125, x.high.cutoff = 3, y.cutoff = 0.8)
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


PrintPCA(object = EXP, pcs.print = 1:5, genes.print = 5, use.full = FALSE)


Stem_gene=c('Prom1','Nes','Egfr','Cd15','Slc1a3','Sox2','Fabp7','Nr2e1','Id3','Clu','Sox9','Vcam1','Slc1a2','Id2','Sox11','Apoe','Tbr2','Ntsr2')

EXP = ScaleData(object = EXP,vars.to.regress = c("percent.mito", "nUMI"), genes.use = EXP@var.genes)
EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )

#EXP = ScaleData(object = EXP,vars.to.regress = c("percent.mito", "nUMI"), genes.use = Stem_gene)
#EXP <- RunPCA(object = EXP, pc.genes = Stem_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )



PCUSE=1:10
PCElbowPlot(object = EXP,num.pc=PCNUM)


PrintPCA(object = EXP, pcs.print = 1:20)



EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE, check_duplicates = FALSE)
#EXP_cluster <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE)


FeaturePlot(object = EXP, features.plot = c('Sox2','Egfr','Id3','Olig2','Zic1','Apoe','Gfap','Slc1a3','Nes'), cols.use = c("grey", "red"), reduction.use = "tsne")
FeaturePlot(object = EXP, features.plot = c('Sox2'),cols.use = c("grey", "red"), reduction.use = "tsne")


#FeaturePlot(object = EXP, features.plot = c('Cd133'),cols.use = c("grey", "red"), reduction.use = "tsne")

save(EXP, file = "EXP.Robj")
save(EXP_cluster, file = "EXP_cluster.Robj")


