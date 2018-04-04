library(Seurat)
library(dplyr)
library(Matrix)
PCNUM=20
PCUSE=1:10
RES=1.2


exp_data=read.table('COMBINED_NORMED.RANK1000.txt',header=T,row.names=1)
EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=200)
EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)
pdf('./VarGene.pdf')
EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.8)
dev.off()
length(x=EXP@var.genes)
EXP = ScaleData(object = EXP,genes.use = EXP@var.genes)
EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE)
#save(EXP, file = "Seurat_EXP.Robj")

EXP_cluster <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE)
#save(EXP_cluster, file = "Seurat_EXP_cluster.Robj")





S10X=which(sapply(colnames(EXP@data), grepl, pattern='S10X'))
GBM=which(sapply(colnames(EXP@data), grepl, pattern='GBM'))

TAG=rep(0,length(colnames(EXP@data)))
TAG[S10X]=1
EXP <- AddMetaData(object = EXP, metadata = TAG, col.name = "TAG")
EXP@meta.data$TAG=as.numeric(TAG)

pdf('BATCH.pdf',width=10,height=10)
FeaturePlot(object = EXP, features.plot = c("TAG"), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()


