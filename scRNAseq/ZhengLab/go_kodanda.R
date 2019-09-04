
setwd('/Volumes/Feng/Zhenglab/Combine')


setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
.set_python('C:/Users/cchmc/Anaconda3/python.exe')

DATA=readRDS('DATA2.RDS') # Kodanda
BATCH=readRDS('BATCH2.RDS')



mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=3, GN=5000, SEED=1, COMBAT=TRUE)
saveRDS(mybeer,file='mybeer2.RDS')




#PCUSE=mybeer$select
PCUSE=.selectUSE(mybeer, CUTR=0.8, CUTL=0.8, RR=0.5, RL=0.5)

COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


pbmc <- mybeer$seurat  
#PCUSE=mybeer$select
PCUSE=.selectUSE(mybeer, CUTR=0.8, CUTL=0.8, RR=0.5, RL=0.5)

#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)

saveRDS(pbmc,file='pbmc2.RDS')






LABEL=readRDS(file='LABEL.RDS')
pbmc@meta.data$celltype=rep(NA,ncol(pbmc))
pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='NATURE')]=LABEL
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)



#######
VEC=pbmc@reductions$umap@cell.embeddings
set.seed(321)
N=100
K=kmeans(VEC,centers=N)
pbmc@meta.data$kclust=K$cluster   



pdf('KDD_CLUST.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='kclust', pt.size=0.1,label=T)+NoLegend()
dev.off()






#'Defa22','Defa24',



pdf('KDD_GENE.pdf',width=10,height=10)
markers.to.plot <- c('Lyz1','Defa17','Ang4','Mmp7')
FeaturePlot(pbmc, features=markers.to.plot)

markers.to.plot <- c('Ptprc','Cd44')
FeaturePlot(pbmc, features=markers.to.plot)

dev.off()

markers.to.plot <- c('Ptprc')


DotPlot(pbmc, features = rev(markers.to.plot), xlab='',ylab='',cols = c("blue", "red"), group.by='kclust',dot.scale = 8) + RotatedAxis()


Idents(pbmc)=pbmc@meta.data$kclust

cluster78.markers <- FindMarkers(pbmc, ident.1 = 78, min.pct = 0.5,only.pos = TRUE)
head(cluster78.markers, n = 5)

cluster95.markers <- FindMarkers(pbmc, ident.1 = 95, min.pct = 0.5,only.pos = TRUE)
head(cluster95.markers, n = 5)





pdf('KDD_GENE.pdf',width=10,height=5)
markers.to.plot <- c('Lyz1','Defa17')
FeaturePlot(pbmc, features=markers.to.plot)

markers.to.plot <- c('Ptprc','Cd44')
FeaturePlot(pbmc, features=markers.to.plot)


markers.to.plot <- c('Slc12a2','Lgr5','Ascl2')
FeaturePlot(pbmc, features=markers.to.plot)

markers.to.plot <- c('Pcna','Cdk4','Clu')
FeaturePlot(pbmc, features=markers.to.plot)

markers.to.plot <- c('Alpi','Fabp1','Apoa1')
FeaturePlot(pbmc, features=markers.to.plot)


dev.off()



Immune=c(31,86,74,38,30,94)
Tuft=c(57)
Goblet=c(77,34,80,49)
Endocrine=c(95,78)
Stem=c(75,58,66)
Paneth=c(72,61,99)
TA=c(79,67,10,60,42,47,62,69,85,56,63,96,82,47,18)

TAG_CLUST=rep('Enterocyte',ncol(pbmc))
TAG_CLUST[which(pbmc@meta.data$kclust %in% Immune)]='Immune'
TAG_CLUST[which(pbmc@meta.data$kclust %in% Tuft)]='Tuft'
TAG_CLUST[which(pbmc@meta.data$kclust %in% Goblet)]='Goblet'
TAG_CLUST[which(pbmc@meta.data$kclust %in% Endocrine)]='Endocrine'
TAG_CLUST[which(pbmc@meta.data$kclust %in% Stem)]='Stem'
TAG_CLUST[which(pbmc@meta.data$kclust %in% Paneth)]='Paneth'
TAG_CLUST[which(pbmc@meta.data$kclust %in% TA)]='TA'



pbmc@meta.data$level1=TAG_CLUST



pdf('KDD_UMAP.pdf',width=6,height=6)
DimPlot(pbmc, reduction.use='umap', group.by='level1', pt.size=0.1,label=T)+NoLegend()
DimPlot(pbmc, cells = colnames(pbmc)[which(pbmc@meta.data$batch=='AGE')],reduction.use='umap', group.by='level1', pt.size=0.1,label=T)+NoLegend()
DimPlot(pbmc, cells = colnames(pbmc)[which(pbmc@meta.data$batch=='YOUNG')],reduction.use='umap', group.by='level1', pt.size=0.1,label=T)+NoLegend()
DimPlot(pbmc, cells = colnames(pbmc)[which(pbmc@meta.data$batch %in% c('AGE','YOUNG'))],reduction.use='umap', group.by='batch', pt.size=0.1,label=T)+NoLegend()
dev.off()




TB=table( pbmc@meta.data$level1, pbmc@meta.data$batch)

TB
.writeTable(TB, PATH='AGEYOUNG_TABLE.txt')


DimPlot(pbmc, reduction.use='umap', group.by='level1', pt.size=0.1,label=T)+NoLegend()
DimPlot(pbmc, cells = colnames(pbmc)[which(pbmc@meta.data$batch=='AGE')],reduction.use='umap', group.by='level1', pt.size=0.1,label=T)+NoLegend()
DimPlot(pbmc, cells = colnames(pbmc)[which(pbmc@meta.data$batch=='YOUNG')],reduction.use='umap', group.by='level1', pt.size=0.1,label=T)+NoLegend()
DimPlot(pbmc, cells = colnames(pbmc)[which(pbmc@meta.data$batch %in% c('AGE','YOUNG'))],reduction.use='umap', group.by='batch', pt.size=0.1,label=T)+NoLegend()



saveRDS(pbmc@meta.data,file='KDD.META.level1.RDS')




############################
USE.CELL=which(pbmc@meta.data$level1=='TA' & pbmc@meta.data$batch!='NATURE')
TAG=pbmc@meta.data$batch
EXP=as.matrix(pbmc@assays$RNA@data[,USE.CELL])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(TAG[USE.CELL]))


OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'./GSEA_KDD/TA_EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'./GSEA_KDD/TA_PT.cls',sep=' ',quote=F,row.names=F,col.names=F )

dim(EXP)




############################
USE.CELL=which(pbmc@meta.data$level1=='Stem' & pbmc@meta.data$batch!='NATURE')
TAG=pbmc@meta.data$batch
EXP=as.matrix(pbmc@assays$RNA@data[,USE.CELL])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(TAG[USE.CELL]))


OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'./GSEA_KDD/Stem_EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'./GSEA_KDD/Stem_PT.cls',sep=' ',quote=F,row.names=F,col.names=F )

dim(EXP)




############################
USE.CELL=which(pbmc@meta.data$level1=='Paneth' & pbmc@meta.data$batch!='NATURE')
TAG=pbmc@meta.data$batch
EXP=as.matrix(pbmc@assays$RNA@data[,USE.CELL])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(TAG[USE.CELL]))


OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'./GSEA_KDD/Paneth_EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'./GSEA_KDD/Paneth_PT.cls',sep=' ',quote=F,row.names=F,col.names=F )

dim(EXP)







