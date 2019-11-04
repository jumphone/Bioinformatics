setwd('F:/HFZ')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

pbmc=readRDS('pbmc.enh.RDS')


VEC=pbmc@reductions$umap@cell.embeddings

# Here, we use K-means to do the clustering
N=100
set.seed(123)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)+NoLegend()




