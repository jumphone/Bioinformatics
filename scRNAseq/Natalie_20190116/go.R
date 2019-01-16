
library('Seurat')
source('scRef.R')

load('Seurat_EXP_cluster.Robj')
com_tag=read.table('2COM.txt',header=T,sep='\t')
inj_tag=read.table('2injury.txt',header=T,sep='\t')
cls_tag=as.numeric(as.character(EXP_cluster@ident))
used_cluster=c(2,9,14,17,19,20,22,23)


pbmc_old=EXP_cluster
used_index=which(cls_tag %in% used_cluster)



pbmc_data=as.matrix(pbmc_old@data[,used_index])
pbmc= CreateSeuratObject(raw.data = pbmc_data, min.cells = 0, min.genes = 0, project = "10X_PBMC")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)




