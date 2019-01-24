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

pdf('CellType.pdf',width=20,height=18)
TSNEPlot(pbmc,group.by='tag',do.label=TRUE)
dev.off()

ref_vec=pbmc@dr$tsne@cell.embeddings
ref_tag=cbind(names(pbmc@ident), as.character(pbmc@meta.data$tag))    
exp_ref_mat=as.matrix(pbmc@raw.data)
LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )  

########################################
scdata=readRDS('R4_umap.RDS')

source('scRef.R')
exp_sc_mat=as.matrix(scdata@raw.data)
  
sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2

out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=300, 
        random_seed=123, CPU=4, print_step=10)

XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))

plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70')
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',)



