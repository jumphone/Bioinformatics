


setwd('/users/zha8dh/tianlab/HFZ')
source('./BEER.R')


pbmc=readRDS('./pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings

# Here, we use K-means to do the clustering
N=50
set.seed(1)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)

saveRDS(pbmc@meta.data, file='META.RDS')


pdf('~/Downloads/HFZ5.CLUST.CHECK.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)
dev.off()


pbmc.markers=readRDS(file='pbmc.markers.RDS')

#EXP_CLUST=.generate_mean(pbmc@assays$RNA@data,pbmc@meta.data$clust)
#######################################
N=300
set.seed(1)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)

saveRDS(pbmc@meta.data, file='META_300.RDS')

pdf('~/Downloads/HFZ5.CLUST.300.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)
dev.off()

EXP_CLUST=.generate_mean(pbmc@assays$RNA@data,pbmc@meta.data$clust)
VEC_CLUST=.generate_mean(t(pbmc@reductions$umap@cell.embeddings), pbmc@meta.data$clust)

saveRDS(EXP_CLUST, file='./EXP_CLUST.RDS')
saveRDS(VEC_CLUST, file='./VEC_CLUST.RDS')


saveRDS(EXP_CLUST, file='~/Downloads/EXP_CLUST.RDS')
saveRDS(VEC_CLUST, file='~/Downloads/VEC_CLUST.RDS')





#####################
