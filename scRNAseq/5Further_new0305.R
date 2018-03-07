library(Seurat)
library(dplyr)
library(Matrix)
load("./images/Seurat_EXP_TSNE.Robj")
load("./images/Seurat_EXP_cluster.Robj")
GFP=read.table('../data_new/GFP/GFP_combined_name.txt',row.names=1)
#GFP=read.table('../data_new/GFP/GFP_combined_nom_name.txt',row.names=1)
#GFP[which(GFP[,1]>0),1]=1

GFP[,1]=log(GFP[,1]+1,2)

colnames(GFP)='GFP'
GFP=data.frame(GFP)

EXP_GFP=AddMetaData(object=EXP,metadata=GFP,col.name='GFP')
EXP_cluster_GFP=AddMetaData(object=EXP_cluster,metadata=GFP,col.name='GFP')


pdf('images/GFP/ALL.pdf',width=10,height=10)
FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne")

dev.off()

S1.5MD=colnames(EXP@data)[which(EXP@ident=='S1.5MD')]
S1.5MP=colnames(EXP@data)[which(EXP@ident=='S1.5MP')]
S1.5MSN=colnames(EXP@data)[which(EXP@ident=='S1.5MSN')]
S4MD=colnames(EXP@data)[which(EXP@ident=='S4MD')]
S4MP=colnames(EXP@data)[which(EXP@ident=='S4MP')]
S4MSN=colnames(EXP@data)[which(EXP@ident=='S4MSN')]


library('ggplot2')
pdf('images/GFP/Samples.pdf',width=10,height=10)
FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne",do.return=T, cells.use =S1.5MD)$GFP +  labs(title = 'S1.5MD')
FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne",do.return=T, cells.use =S1.5MP)$GFP +  labs(title = 'S1.5MP')
FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne",do.return=T, cells.use =S1.5MSN)$GFP  + labs(title = 'S1.5MSN')
FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne",do.return=T, cells.use =S4MD)$GFP  + labs(title = 'S4MD')
FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne",do.return=T, cells.use =S4MP)$GFP  + labs(title = 'S4MP')
FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne",do.return=T, cells.use =S4MSN)$GFP  + labs(title = 'S4MSN')
dev.off()

pdf('images/GFP/Clusters.pdf',width=10,height=10)
P=c()
i=0
while(i <=30){
CELL=colnames(EXP_cluster@data)[which(EXP_cluster@ident==i)]
p=FeaturePlot(object = EXP_cluster_GFP, features.plot = c("GFP"), cols.use = c("grey90", "red"), reduction.use = "tsne",do.return=T, cells.use =CELL)$GFP + ggtitle(as.character(i))
ggsave(paste0("images/GFP/Cluster/",as.character(i),".png"), dpi=300)
print(i)
i=i+1}
dev.off()


pdf('images/GFP/EXP.pdf',width=30,height=30)
FeaturePlot(object = EXP_cluster, features.plot = c("Acta2","Itgb1","S100a4","Cav1","Fap","Tnfsf4","Jam2"), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Foxp3","Cd4","Ccl5","Car2","Ptprc"), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Cxcr4","Cxcl12","Cxcr3","Cxcl10","Foxr2","Tnfrsf4","Pdcd1lg2","Pdgfrl","Cd274"), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()


#cell number in an excel file for GFP+ cells in each cluster for all 6 samples
i=0
while(i<=30){
nonGFP=table(EXP_GFP@ident[which(EXP_cluster@ident==i & EXP_GFP@meta.data$GFP==0)])
GFP=table(EXP_GFP@ident[which(EXP_cluster@ident==i & EXP_GFP@meta.data$GFP>0)])
c=cbind(nonGFP,GFP)
print(c)
write.table(c,file=paste0('images/GFP/CNUM/',as.character(i),'.csv'),row.names=T,col.names=T,quote=F,sep='\t')
i=i+1
}
i=0
nonGFP=table(EXP_GFP@ident[which(EXP_cluster@ident==i & EXP_GFP@meta.data$GFP==0)])
GFP=table(EXP_GFP@ident[which(EXP_cluster@ident==i & EXP_GFP@meta.data$GFP>0)])
c=cbind(nonGFP,GFP)

i=1
while(i<=30){
nonGFP=table(EXP_GFP@ident[which(EXP_cluster@ident==i & EXP_GFP@meta.data$GFP==0)])
GFP=table(EXP_GFP@ident[which(EXP_cluster@ident==i & EXP_GFP@meta.data$GFP>0)])
c_new=cbind(nonGFP,GFP)
c=cbind(c,c_new)
#print(c)
#write.table(c,file=paste0('images/GFP/CNUM/',as.character(i),'.tsv'),row.names=T,col.names=T,quote=F,sep='\t')
i=i+1
}

write.table(c,file=paste0('images/GFP/CNUM_GFP.tsv'),row.names=T,col.names=T,quote=F,sep='\t')



S1.5GFP=which(EXP_GFP@meta.data$GFP > 0 & EXP_GFP@ident %in% c('S1.5MD','S1.5MP','S1.5MSN'))
S4GFP=which(EXP_GFP@meta.data$GFP > 0 & EXP_GFP@ident %in% c('S4MD','S4MP','S4MSN'))


write.table(as.matrix(EXP_GFP@data)[,S1.5GFP],file='images/GFP/S1.5GFP.tsv',row.names=T,col.names=T,quote=F,sep='\t')
write.table(as.matrix(EXP_GFP@data)[,S4GFP],file='images/GFP/S4GFP.tsv',row.names=T,col.names=T,quote=F,sep='\t')



