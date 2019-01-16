
library('Seurat')
source('scRef.R')

load('Seurat_EXP_cluster.Robj')

com_tag=read.table('2COM.txt',header=T,sep='\t')
inj_tag=read.table('2injury.txt',header=T,sep='\t')
all_tag=read.table('cMCA.txt',header=T,sep='\t')
cls_tag=as.numeric(as.character(EXP_cluster@ident))

used_cluster=c(2,9,14,17,19,23)
pbmc_old=EXP_cluster
used_index=which(cls_tag %in% used_cluster & (!com_tag[,2] %in% c('Macrophage_Lyz2.high.Brain.','B.cell_Igkc.high.Bone.Marrow.','T.cell_Ms4a4b.high.Bone.Marrow.')))


pbmc_data=as.matrix(pbmc_old@data[,used_index])
pbmc= CreateSeuratObject(raw.data = pbmc_data, min.cells = 0, min.genes = 0, project = "10X_PBMC")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
#VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
tmp_data=pbmc@data
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc@data=tmp_data
pbmc <- FindVariableGenes(object = pbmc,do.plot=FALSE, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
allgene=rownames(pbmc@data)
#4901
#############
pbmc@meta.data$cls=cls_tag[used_index]
pbmc@meta.data$com=com_tag[used_index,2]
pbmc@meta.data$inj=inj_tag[used_index,2]
pbmc@meta.data$all=as.character(all_tag[used_index,2])
#############
pbmc@meta.data$allm=pbmc@meta.data$all
lname=names(table(pbmc@meta.data$allm))[which(table(pbmc@meta.data$allm)<50)]
pbmc@meta.data$allm[which(pbmc@meta.data$allm %in% lname)]='NA'
#############
pbmc@meta.data$age=pbmc@meta.data$orig.ident
############
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
PCNUM=150
#pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = allgene, do.print = FALSE, pcs.print = 1:5, genes.print = 5)
#PCElbowPlot(object = pbmc)
####

####
PCUSE=1:35

#pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)

#DimPlot(pbmc, reduction.use = "tsne",group.by='cls',pt.size=0.5,do.label=TRUE)
DimPlot(pbmc, reduction.use = "tsne",group.by='age',pt.size=0.5,do.label=TRUE)
#DimPlot(pbmc, reduction.use = "tsne",group.by='com',pt.size=0.5,do.label=TRUE)
#DimPlot(pbmc, reduction.use = "tsne",group.by='inj',pt.size=0.5,do.label=TRUE)

#pbmc <- RunUMAP(object = pbmc, dims.use = PCUSE)
DimPlot(pbmc, reduction.use = "umap",group.by='cls',pt.size=0.5,do.label=TRUE)
DimPlot(pbmc, reduction.use = "umap",group.by='age',pt.size=0.5,do.label=TRUE)
DimPlot(pbmc, reduction.use = "umap",group.by='com',pt.size=0.5,do.label=TRUE)
DimPlot(pbmc, reduction.use = "umap",group.by='inj',pt.size=0.5,do.label=TRUE)
DimPlot(pbmc, reduction.use = "umap",group.by='allm',pt.size=0.5,do.label=TRUE)
