library(Seurat)
library(dplyr)
library(Matrix)

load('SSN_EXP_cluster.Robj')
SSN_CLUSTER=EXP_cluster@ident

#load('../../images/Seurat_EXP_cluster.Robj')
#OLD_CLUSTER=EXP_cluster@ident
#EXP_cluster@ident=SSN_CLUSTER
exp_scale = t(apply(EXP_cluster@data,1,scale,scale=T,center=T))

Compare = function(CA,CB){
A=c()
B=c()
Z=c()
P=c()
i=1
while(i <=length(exp_scale[,1])){
a=t(exp_scale[i,which(EXP_cluster@ident %in% as.factor(CA))])
b=t(exp_scale[i,which(EXP_cluster@ident %in% as.factor(CB))])
#pvalue=wilcox.test(a,b) #min(pnorm(z),1-pnorm(z))*2
pvalue=t.test(a,b) #min(pnorm(z),1-pnorm(z))*2
P=c(P,pvalue$p.value)
#print(i)
i=i+1
}
return(P)
}

i=0
while(i<=length(table(EXP_cluster@ident))){
print(i)
#TSNEPlot(object = EXP_cluster,colors.use=c(rep('grey',i),'red',rep('grey',30-i)))
#cluster.markers <- FindMarkers(object = EXP_cluster, ident.1 = i, min.pct = 0.25)
#cluser_top = head(x = cluster.markers, n = 1000)
P=Compare(c(i),c(0:length(table(EXP_cluster@ident))))
AP=p.adjust(P,method='bonferroni')
cluser_top= as.matrix(cbind(row.names(EXP_cluster@data),P,AP))
write.table(file=paste0('Cluster_Marker_Pair/Cluster_',as.character(i),'_marker.tsv'),cluser_top,sep='\t',quote=F,row.names=F,col.names=F)
i=i+1}
