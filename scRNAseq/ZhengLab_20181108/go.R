#/ZhengLab/10X_20181030/OUT
library(Seurat)

###### Combine all files ###################
WT291_1.data <- Read10X(data.dir = "./10X_291_WT_20180925_3mm_1_out/filtered_gene_bc_matrices/mm10/")
WT291_2.data <- Read10X(data.dir = "./10X_291_WT_20180925_3mm_2_out/filtered_gene_bc_matrices/mm10/")
WT291_3.data <- Read10X(data.dir = "./10X_291_WT_20180925_3mm_3_out/filtered_gene_bc_matrices/mm10/")
WT291_4.data <- Read10X(data.dir = "./10X_291_WT_20180925_3mm_4_out/filtered_gene_bc_matrices/mm10/")

KO294_1.data <- Read10X(data.dir = "./10X_294_KO_20180925_3mm_1_out/filtered_gene_bc_matrices/mm10/")
KO294_2.data <- Read10X(data.dir = "./10X_294_KO_20180925_3mm_2_out/filtered_gene_bc_matrices/mm10/")
KO294_3.data <- Read10X(data.dir = "./10X_294_KO_20180925_3mm_3_out/filtered_gene_bc_matrices/mm10/")
KO294_4.data <- Read10X(data.dir = "./10X_294_KO_20180925_3mm_4_out/filtered_gene_bc_matrices/mm10/")


WT291_1=CreateSeuratObject(raw.data = WT291_1.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")
WT291_2=CreateSeuratObject(raw.data = WT291_2.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")
WT291_3=CreateSeuratObject(raw.data = WT291_3.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")
WT291_4=CreateSeuratObject(raw.data = WT291_4.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")


KO294_1=CreateSeuratObject(raw.data = KO294_1.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")
KO294_2=CreateSeuratObject(raw.data = KO294_2.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")
KO294_3=CreateSeuratObject(raw.data = KO294_3.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")
KO294_4=CreateSeuratObject(raw.data = KO294_4.data, min.cells = 0, min.genes = 0,  project = "ZhengLab")

WT12 = MergeSeurat(WT291_1, WT291_2, project = NULL, min.cells = 0, min.genes = 0, is.expr = 0, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE, names.field = 1, names.delim = "_", add.cell.id1 = 'WT1', add.cell.id2 = 'WT2')
WT34 = MergeSeurat(WT291_3, WT291_4, project = NULL, min.cells = 0, min.genes = 0, is.expr = 0, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE, names.field = 1, names.delim = "_", add.cell.id1 = 'WT3', add.cell.id2 = 'WT4')
WT = MergeSeurat(WT12, WT34, project = NULL, min.cells = 0, min.genes = 0, is.expr = 0, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE, names.field = 1,names.delim = "_", add.cell.id1 = '', add.cell.id2 = '')

KO12 = MergeSeurat(KO294_1, KO294_2, project = NULL, min.cells = 0, min.genes = 0, is.expr = 0, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE, names.field = 1, names.delim = "_", add.cell.id1 = 'KO1', add.cell.id2 = 'KO2')
KO34 = MergeSeurat(KO294_3, KO294_4, project = NULL, min.cells = 0, min.genes = 0, is.expr = 0, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE, names.field = 1, names.delim = "_", add.cell.id1 = 'KO3', add.cell.id2 = 'KO4')
KO = MergeSeurat(KO12, KO34, project = NULL, min.cells = 0, min.genes = 0, is.expr = 0, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE, names.field = 1,names.delim = "_", add.cell.id1 = '', add.cell.id2 = '')

EXP = MergeSeurat(WT, KO, project = NULL, min.cells = 0, min.genes = 0, is.expr = 0, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE, names.field = 1,names.delim = "_", add.cell.id1 = 'WT', add.cell.id2 = 'KO')

save(EXP,file='EXP.RData')

##################################################################

library(Seurat)
load('EXP.RData')
pbmc=EXP
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(1000, -Inf), high.thresholds = c(2000, 0.2))

table(pbmc@ident)

