library(Seurat)
library(dplyr)
library(Matrix)

load('./images/Seurat_EXP_cluster.Robj')

pdf('./images/Seurat_TSNE_clusters.pdf',width=30,height=15)
#TSNEPlot(object = EXP_cluster,colors.use=c(rep('grey',11),'red',rep('grey',19)))
i=0
while(i<=30){
print(i)
TSNEPlot(object = EXP_cluster,colors.use=c(rep('grey',i),'red',rep('grey',30-i)))
cluster.markers <- FindMarkers(object = EXP_cluster, ident.1 = i, min.pct = 0.25)
cluser_top = head(x = cluster.markers, n = 1000)
write.table(file=paste0('./images/Cluster_',as.character(i),'_marker.tsv'),cluser_top,sep='\t',quote=F)
i=i+1}
dev.off()

