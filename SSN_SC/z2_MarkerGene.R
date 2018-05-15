library(Seurat)
library(dplyr)
library(Matrix)

load('SSN_EXP_cluster.Robj')
SSN_CLUSTER=EXP_cluster@ident
SSN_CELL=colnames(EXP_cluster@data)

load('../../images/Seurat_EXP_cluster.Robj')
OLD_CLUSTER=EXP_cluster@ident
OLD_CELL=colnames(EXP_cluster@data)
TMP=as.character(OLD_CLUSTER)
SSN_CLUSTER=as.character(SSN_CLUSTER)
TMP[which(OLD_CELL %in% SSN_CELL)]=SSN_CLUSTER
NAF=as.character(max(as.numeric(SSN_CLUSTER))+1)
TMP[which(! OLD_CELL %in% SSN_CELL)]=NAF

TMP=as.factor(TMP)
names(TMP)=names(EXP_cluster@ident)
EXP_cluster@ident=TMP




i=0
while(i<=length(table(EXP_cluster@ident))){
print(i)
#TSNEPlot(object = EXP_cluster,colors.use=c(rep('grey',i),'red',rep('grey',30-i)))
cluster.markers <- FindMarkers(object = EXP_cluster, ident.1 = i, min.pct = 0.25)
cluser_top = head(x = cluster.markers, n = 1000)
write.table(file=paste0('Cluster_Marker_Gene/Cluster_',as.character(i),'_marker.tsv'),cluser_top,sep='\t',quote=F)
i=i+1}
