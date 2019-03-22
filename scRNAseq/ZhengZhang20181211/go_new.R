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
    low.thresholds = c(200, -Inf), high.thresholds = c(7000, 1))
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
pdf('PCA.pdf')
PCElbowPlot(object = pbmc)
dev.off()

PCUSE=1:20
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
TSNEPlot(object = pbmc,pt.size=0.5)

pdf('TSNE.pdf')
TSNEPlot(object = pbmc,pt.size=0.5)
dev.off()


#saveRDS(pbmc,file='WT.RDS')
#######################


#pbmc=readRDS('WT.RDS')

##

library(Seurat)
#source('scRef.R')


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]



exp_ref_mat=read.table('MCA_INTEST.txt',header=T,row.names=1,sep='\t')
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
pbmc@meta.data$mca=out$tag2[,2]

pdf('CellType.pdf',width=14,height=10)
TSNEPlot(object = pbmc, do.label=T, group.by ='mca', pt.size = 0.5)
dev.off()






saveRDS(pbmc,file='WT.RDS')

########################################################################

pbmc=readRDS('WT.RDS')

plot(pbmc@meta.data$nGene, pbmc@meta.data$nUMI)
VEC=cbind(pbmc@meta.data$nGene, pbmc@meta.data$nUMI)
K=kmeans(VEC,centers=3)
C=K$cluster

plot(pbmc@meta.data$nGene, pbmc@meta.data$nUMI,col=C,pch=16)
table(C)
table(C)/sum(table(C))
saveRDS(K,'K.RDS')
#0.931182796 0.002977667 0.065839537
############




library(Seurat)
pbmc=readRDS('WT.RDS')
K=readRDS('K.RDS')
C=K$cluster

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

exp_sc_mat=exp_sc_mat[,which(C==1)]



source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
#pbmc.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(raw.data = exp_sc_mat, min.cells = 3, min.genes = 200, 
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
    low.thresholds = c(200, -Inf), high.thresholds = c(7000, 1))
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



#######################

##

library(Seurat)
#source('scRef.R')


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]



exp_ref_mat=read.table('MCA_INTEST.txt',header=T,row.names=1,sep='\t')
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
pbmc@meta.data$mca=out$tag2[,2]

pdf('new_CellType.pdf',width=14,height=10)
TSNEPlot(object = pbmc, do.label=T, group.by ='mca', pt.size = 0.5)
dev.off()



#saveRDS(pbmc,file='new_WT.RDS')
