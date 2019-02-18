library(Seurat)
library(dplyr)

# Load the PBMC dataset
require(data.table)
pbmc.data=data.frame(fread("E6701_exp.txt",header=T), row.names=1)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 100, min.genes = 1000, project = "10X_PBMC")

rm(pbmc.data)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(5000, 0.1))

#V=apply(pbmc@data,1,var)

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F, mean.function = ExpMean, dispersion.function = LogVMR,  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI", "percent.mito"))


PCNUM=50
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, pcs.compute=PCNUM, genes.print = 5)
PCElbowPlot(object = pbmc,num.pc=PCNUM)




PCUSE=1:50
#pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pbmc <- RunUMAP(pbmc, dims.use = PCUSE)
DimPlot(pbmc, reduction.use = "umap",pt.size=0.1)


tag=read.table('E6701_cluster.txt',sep='\t')
tag=tag[which(tag[,1] %in% colnames(pbmc@data)),]

pbmc@meta.data$F=tag[,2]
pbmc@meta.data$C=tag[,3]
pbmc@meta.data$L=tag[,4]

pdf('13Nature.pdf',width=12,height=8)
DimPlot(pbmc, reduction.use = "umap",pt.size=0.1,group.by='F')
DimPlot(pbmc, reduction.use = "umap",pt.size=0.1,group.by='C')
DimPlot(pbmc, reduction.use = "umap",pt.size=0.1,group.by='L')
dev.off()

#saveRDS(pbmc,file='13Nature_10X.RDS')


#########################
#pbmc=readRDS('13Nature_10X.RDS')

source('scRef.R')
exp_ref_mat=read.table('MCA_combined_human.txt',header=T,row.names=1,sep='\t')


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

out=SCREF(exp_sc_mat, exp_ref_mat, min_cell=100, CPU=4,print_step=10)
tag1=out$tag1
tag2=out$tag2

pbmc@meta.data$scref=tag2[,2]
pbmc@meta.data$tmp=tag1[,2]

pdf('placenta_MCA.pdf',width=12,height=10)
DimPlot(pbmc, reduction.use = "umap",group.by='scref', do.label=T, pt.size=3, label.size=5)
dev.off()










