library(Seurat)
library(dplyr)

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
pbmc.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pdf('WT_QC.pdf')
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(4000, 0.5))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 1, y.cutoff = 0.1)
length(x = pbmc@var.genes)

#all.genes=rownames(pbmc@raw.data)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

PCNUM=50
#pbmc <- RunPCA(object = pbmc, pc.genes = all.genes,pcs.compute = PCNUM, do.print = F, pcs.print = 1:5, 
#    genes.print = 5)
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes,pcs.compute = PCNUM, do.print = F, pcs.print = 1:5, 
    genes.print = 5)
pdf('WT_PCA.pdf')
PCElbowPlot(object = pbmc)
dev.off()

PCUSE=1:20
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
TSNEPlot(object = pbmc,pt.size=0.5)

pdf('WT_TSNE.pdf')
TSNEPlot(object = pbmc,pt.size=0.5)
dev.off()


#saveRDS(pbmc,file='WT.RDS')
#######################


#pbmc=readRDS('WT.RDS')

##

library(Seurat)
source('scRef.R')


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]



exp_ref_mat=read.table('MCA_SmallIntestine_ref_mouse.txt',header=T,row.names=1,sep='\t')
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
pbmc@meta.data$MCASmallIntest=out$tag2[,2]

pdf('WT_MCASmallIntest.pdf',width=14,height=10)
TSNEPlot(object = pbmc, do.label=T, group.by ='MCASmallIntest', pt.size = 0.5)
dev.off()




exp_ref_mat=read.table('TM_Large_Intestine_ref_mouse.txt',header=T,row.names=1,sep='\t')
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[3])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
pbmc@meta.data$TMLargeIntest=out$tag2[,2]

pdf('WT_TMLargeIntest.pdf',width=14,height=10)
TSNEPlot(object = pbmc, do.label=T, group.by ='TMLargeIntest', pt.size = 0.5)
dev.off()


saveRDS(pbmc,file='WT.RDS')

########################################################################



library(Seurat)
library(dplyr)

source('scRef.R')
pbmc.data <- Read10X(data.dir = "./KO/mm10/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pdf('KO_QC.pdf')
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    low.thresholds = c(200, -Inf), high.thresholds = c(7000, 0.5))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 1, y.cutoff = 0.1)
length(x = pbmc@var.genes)

#all.genes=rownames(pbmc@raw.data)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

PCNUM=50
#pbmc <- RunPCA(object = pbmc, pc.genes = all.genes,pcs.compute = PCNUM, do.print = F, pcs.print = 1:5, 
#    genes.print = 5)
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes,pcs.compute = PCNUM, do.print = F, pcs.print = 1:5, 
    genes.print = 5)
pdf('KO_PCA.pdf')
PCElbowPlot(object = pbmc)
dev.off()

PCUSE=1:20
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
TSNEPlot(object = pbmc,pt.size=0.5)

pdf('KO_TSNE.pdf')
TSNEPlot(object = pbmc,pt.size=0.5)
dev.off()


saveRDS(pbmc,file='KO.RDS')
#######################


#pbmc=readRDS('WT.RDS')

##

library(Seurat)
source('scRef.R')


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]



exp_ref_mat=read.table('MCA_SmallIntestine_ref_mouse.txt',header=T,row.names=1,sep='\t')
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
pbmc@meta.data$MCASmallIntest=out$tag2[,2]

pdf('KO_MCASmallIntest.pdf',width=14,height=10)
TSNEPlot(object = pbmc, do.label=T, group.by ='MCASmallIntest', pt.size = 1)
dev.off()




exp_ref_mat=read.table('TM_Large_Intestine_ref_mouse.txt',header=T,row.names=1,sep='\t')
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[3])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
pbmc@meta.data$TMLargeIntest=out$tag2[,2]

pdf('KO_TMLargeIntest.pdf',width=14,height=10)
TSNEPlot(object = pbmc, do.label=T, group.by ='TMLargeIntest', pt.size = 0.5)
dev.off()


#saveRDS(pbmc,file='KO.RDS')
