library(Seurat)

PCNUM=20
PCUSE=1:10
RES=0.8

exp_data=read.table('run1642_10000.dge.txt',header=T,row.names=1)
EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=200)

mito.genes <- grep(pattern = "^Mt", x = rownames(x = EXP@data), value = TRUE)
percent.mito <- colSums(EXP@data[mito.genes, ]) / colSums(EXP@data)
EXP <- AddMetaData(object = EXP, metadata = percent.mito, col.name = "percent.mito")


pdf('Seurat_QC.pdf',width=30,height=15)
VlnPlot(object = EXP, features.plot = c("nGene", "nUMI",'percent.mito'), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "nGene")
dev.off()

EXP=FilterCells(object = EXP, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(3000, 0.05))

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

pdf('Seurat_VarGene.pdf')
EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.3)
dev.off()

length(x=EXP@var.genes)

'Prom1' %in% rownames(EXP@data)
'Prom1' %in% EXP@var.genes
'Nes' %in% EXP@var.genes
'Egfr' %in% EXP@var.genes
'Cd15' %in% EXP@var.genes
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
