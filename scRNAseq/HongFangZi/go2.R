


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





