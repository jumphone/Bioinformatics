

library(dplyr)
library(Seurat)

mybeer=readRDS('mybeer.RDS')
PCUSE <- which(mybeer$cor>quantile(mybeer$cor,0.3) )
#npbmc <- mybeer$seurat
#npbmc <- RunUMAP(object = npbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
#DimPlot(npbmc, reduction = "umap")


pbmc=readRDS('ALL.RDS')

STEM=c("Lgr5", "Ascl2", "Slc12a2", "Axin2", "Olfm4", "Axin2")
TA=c("Mki67", "Cdk4", "Mcm5", "Mcm6", "Pcna")
QSTEM=c("Bmi1", "Hopx", "Lrig1", "Tert")



UMAP=pbmc@reductions$umap@cell.embeddings
set.seed(12345)
KM=kmeans(UMAP,centers=30)
KMC=KM$cluster
plot(UMAP,col=KMC,pch=15)
CLUST=KMC

CLUST=as.factor(CLUST)
names(CLUST)=names(pbmc@active.ident)
pbmc@active.ident=CLUST
DimPlot(pbmc,label=T)


pbmc@meta.data$clust=as.character(CLUST)

pdf('EXP_DotPlot.pdf',width=12,height=8)

DimPlot(pbmc,group.by='clust',label=T)

markers.to.plot <- c("Lgr5","Olfm4","Ascl2","Axin2","Mki67",'Mcm5',"Cdk4",'Pcna')
DotPlot(pbmc, features = rev(markers.to.plot), xlab='',ylab='',cols = c("blue", "red"), group.by='clust',dot.scale = 8) + RotatedAxis()


FeaturePlot(pbmc, features=STEM)
FeaturePlot(pbmc, features=TA)
VlnPlot(pbmc,group.by='clust',features=c('Cdk2'))
VlnPlot(pbmc,group.by='clust',features=c('Cdk4'))
VlnPlot(pbmc,group.by='clust',features=c('Top2a'))
VlnPlot(pbmc,group.by='clust',features=c('Mki67'))
dev.off()


pbmc@meta.data$newall=pbmc@meta.data$clust
pbmc@meta.data$newall[which(pbmc@meta.data$clust=='7')]='Stem'
pbmc@meta.data$newall[which(pbmc@meta.data$clust %in% c('14','16'))]='TA'
pbmc@meta.data$newall[which(pbmc@meta.data$clust %in% c('6','15','8'))]='ImmuneCell'
pbmc@meta.data$newall[which(!pbmc@meta.data$newall %in% c('Stem','TA','ImmuneCell'))]='Enterocyte'

DimPlot(pbmc, group.by='newall',label=T)
#source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
source('scRef.R')

USED=which(pbmc@meta.data$newall %in% c('Enterocyte'))
exp_sc=as.matrix(npbmc@assays$RNA@data)
exp_sc=exp_sc[,USED]
exp_ref=read.table('GSE92332_intestin_mouse_ref.txt',header=T,row.names=1,sep='\t')
exp_ref=exp_ref[,which(colnames(exp_ref) %in% c('Enterocyte.Immature.Distal','Enterocyte.Immature.Proximal',
                                                'Enterocyte.Mature.Distal','Enterocyte.Mature.Proximal',
                                               'Enterocyte.Progenitor'))]
exp_ref[,1]=exp_ref[,1]+exp_ref[,2]
exp_ref[,3]=exp_ref[,3]+exp_ref[,4]
exp_ref=exp_ref[,c(1,3,5)]
colnames(exp_ref)=c('Enterocyte.Immature','Enterocyte.Mature','Enterocyte.Progenitor')

OUT=SCREF(exp_sc, exp_ref, CPU=4, min_cell=10,print_step=10)

pbmc@meta.data$newall[USED]=OUT$tag2[,2]
pbmc@meta.data$newall[which(pbmc@meta.data$newall =='Enterocyte.Immature.Distal')]='Enterocyte.Immature'
pbmc@meta.data$newall[which(pbmc@meta.data$newall =='Enterocyte.Mature.Distal')]='Enterocyte.Mature'




pdf('NEWTYPE.pdf',width=12,height=8)
DimPlot(pbmc, group.by='newall',label=T)
VlnPlot(pbmc,group.by='newall',features=c('Cdk2'))
VlnPlot(pbmc,group.by='newall',features=c('Cdk4'))
VlnPlot(pbmc,group.by='newall',features=c('Top2a'))
VlnPlot(pbmc,group.by='newall',features=c('Mki67'))
VlnPlot(pbmc,group.by='newall',features=c('Mki67'))
dev.off()



Idents(pbmc)=pbmc@meta.data$newall
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

save.image('TRY1.RData')
