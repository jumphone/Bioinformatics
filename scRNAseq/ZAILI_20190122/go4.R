library(Seurat)
library(dplyr)

pbmc.data <- readRDS('cerebullum_dev.RDS')
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
      
length(x = pbmc@var.genes)
#4602

pbmc <- ScaleData(object = pbmc, genes.use =pbmc@var.genes, vars.to.regress = c("nUMI", "percent.mito"),num.cores =4,do.par=TRUE)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

meta.data=read.table('cerebellum_cell_metadata.tsv',sep='\t',header=T)

pbmc@dr$tsne@cell.embeddings[,1]=as.numeric(meta.data[,1])
pbmc@dr$tsne@cell.embeddings[,2]=as.numeric(meta.data[,2])
pbmc@meta.data$tag=as.character(meta.data[,15])
pbmc@meta.data$dev=as.character(meta.data[,4])
pbmc@meta.data$clust=meta.data[,14]

saveRDS(pbmc,'cb_seurat.RDS')
########################################

pbmc=readRDS('cb_seurat.RDS')
pbmc@meta.data$clust=paste0('C',pbmc@meta.data$clust)

TAB=table(pbmc@meta.data$clust,pbmc@meta.data$tag)
write.table(TAB,file='Sup1_CellTypeAndClust.txt',sep='\t',quote=F,row.names=T,col.names=T)

pdf('Sup2_CellTypeAndClust.pdf',width=20,height=18)
TSNEPlot(pbmc,group.by='tag',do.label=TRUE)
TSNEPlot(pbmc,group.by='clust',do.label=TRUE)
dev.off()

source('scRef.R')
ref_vec=pbmc@dr$tsne@cell.embeddings
ref_tag=cbind(names(pbmc@ident), as.character(pbmc@meta.data$clust))    
exp_ref_mat=as.matrix(pbmc@raw.data)
rownames(exp_ref_mat)=toupper(rownames(exp_ref_mat))

#####remove dup gene name###
tmp_tab=table(rownames(exp_ref_mat))
tmp_gene=names(which(tmp_tab>1))
exp_ref_mat=exp_ref_mat[which(!rownames(exp_ref_mat) %in% tmp_gene),]
########################################

LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )  
write.table(LocalRef,'Sup3_ClustRef.txt',sep='\t',quote=F,row.names=T,col.names=T)
saveRDS(exp_ref_mat,'./Sup4_projection/exp_ref_mat.RDS')
saveRDS(ref_tag,'./Sup4_projection/ref_tag.RDS')
saveRDS(ref_vec,'./Sup4_projection/ref_vec.RDS')


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

