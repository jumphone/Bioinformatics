
library(Seurat)
source('scRef.R')

tag_data=read.table('GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt',sep='\t',row.names=1,header=T)


pbmc.data=read.table('GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt',sep='\t',row.names=1,header=T)
#saveRDS(pbmc.data,'ALL.RDS')

MS=pbmc.data[,which(tag_data[,3]=='MS')]
CT=pbmc.data[,which(tag_data[,3]=='Ctrl')]

saveRDS(MS,'MS.RDS')
saveRDS(CT,'CT.RDS')

CPU=6
PCNUM=50
PCUSE=1:PCNUM

##########
MSD = CreateSeuratObject(raw.data = MS, min.cells = 0, min.genes = 0, project = "MS") 
MSD <- NormalizeData(object = MSD, normalization.method = "LogNormalize", scale.factor = 10000)
MSD <- ScaleData(object = MSD, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
MSD <- RunPCA(object = MSD, pcs.compute=PCNUM, pc.genes = rownames(MSD@data), do.print = FALSE)
MSD <- RunTSNE(object = MSD, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
MSX=MSD@dr$tsne@cell.embeddings
saveRDS(MSX,file='MSX.RDS')
##########
CTD = CreateSeuratObject(raw.data = CT, min.cells = 0, min.genes = 0, project = "CT") 
CTD <- NormalizeData(object = CTD, normalization.method = "LogNormalize", scale.factor = 10000)
CTD <- ScaleData(object = CTD, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
CTD <- RunPCA(object = CTD, pcs.compute=PCNUM, pc.genes = rownames(MSD@data), do.print = FALSE)
CTD <- RunTSNE(object = CTD, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
CTX=CTD@dr$tsne@cell.embeddings
saveRDS(CTX,file='CTX.RDS')
##########
##########
MSX=readRDS(MSX)
CTX=readRDS(CTX)




saveRDS(pbmc,file='CCA.RDS')




library(kBET)
batch.estimate <- kBET(as.matrix(pbmc@data), pbmc@meta.data$paper.condition)


saveRDS(pbmc,file='CCA.RDS')

DimPlot(pbmc, group.by='orig.ident', reduction.use='tsne')
DimPlot(pbmc, group.by='paper.type', reduction.use='tsne')
DimPlot(pbmc, group.by='paper.cluster', reduction.use='tsne')
DimPlot(pbmc, group.by='paper.condition', reduction.use='tsne')





#pbmc <- RunUMAP(object = pbmc, dims.use = PCUSE)
#DimPlot(pbmc, group.by='orig.ident', reduction.use='umap')
#pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = PCUSE, 
#    resolution = 0.6, print.output = 0, save.SNN = TRUE)


