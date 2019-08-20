


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
N=50
set.seed(1)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)

#saveRDS(pbmc@meta.data, file='META_300.RDS')

#pdf('~/Downloads/HFZ5.CLUST.300.pdf',width=20,height=20)
#DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE) + NoLegend()
#dev.off()

EXP_CLUST=.generate_mean(pbmc@assays$RNA@data,pbmc@meta.data$clust)
VEC_CLUST=.generate_mean(t(pbmc@reductions$umap@cell.embeddings), pbmc@meta.data$clust)

saveRDS(EXP_CLUST, file='./EXP_CLUST.RDS')
saveRDS(VEC_CLUST, file='./VEC_CLUST.RDS')


saveRDS(EXP_CLUST, file='~/Downloads/EXP_CLUST.RDS')
saveRDS(VEC_CLUST, file='~/Downloads/VEC_CLUST.RDS')


##########################
pbmc.markers=readRDS(file='pbmc.markers.RDS')
.writeTable(pbmc.markers,PATH='~/Downloads/Markers.txt')

library(dplyr)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10=as.matrix(top10)
rownames(top10)=c(1:nrow(top10))
.writeTable(top10,PATH='~/Downloads/top10.txt')


pbmc=readRDS('./pbmc.RDS')
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pbmc@meta.data=readRDS(file='META.RDS')
Idents(pbmc)=pbmc@meta.data$clust

pdf('~/Downloads/HFZ6.HEAT.pdf',width=30,height=30)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()
#####################














