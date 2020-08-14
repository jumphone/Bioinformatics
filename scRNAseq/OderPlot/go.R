pseudoExp<-function (pbmc,  gene, ps, clust){
  #ind=which(rownames(pbmc@assays$RNA@data)==gene)
  ind=which(pbmc$cell.type %in% clust)
  ind=rownames(pbmc@meta.data)[ind]
  used_exprs=pbmc@assays$RNA@data[gene,ind]
  
  ind2=which(rownames(ps)%in% ind)
  ps2=log(ps+1)
  used_ps2=ps2[ind2]
  used_order=order(used_ps2)
  all=cbind(used_ps2[used_order], used_exprs[used_order])
  
  #all[,2]=sample(all[,2])
  all[,2]=smooth(all[,2])
  plot(all[,1],all[,2], xlab="Log Pseudotime",type='both',
       pch=19,ylab="Normalized Expression")  
}
    
setwd('F:/OderPlot')
#library('Seurat')
    
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
 pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")   
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)    
   pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

 all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)   
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    
    
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
    pbmc <- RunUMAP(pbmc, dims = 1:10)
    
    
    saveRDS(pbmc,file='pbmc.rds')
########################

ps=rank(pbmc@assays$RNA@data[which(rownames(pbmc)=='TP53'),])
ps=as.matrix(ps,ncol=1)
rownames(ps)=rownames(pbmc@meta.data)
clust=2
gene='TP53'
pbmc$cell.type=pbmc@meta.data$seurat_clusters
pseudoExp(pbmc,  gene, ps, clust)    
    
