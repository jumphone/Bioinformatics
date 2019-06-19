

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

VlnPlot(pbmc,group.by='clust',features=c('Mcm5'))
VlnPlot(pbmc,group.by='clust',features=c('Mcm6'))
VlnPlot(pbmc,group.by='clust',features=c('Mki67'))
VlnPlot(pbmc,group.by='clust',features=c('Bmi1'))
dev.off()


pbmc@meta.data$newall=pbmc@meta.data$clust
pbmc@meta.data$newall[which(pbmc@meta.data$clust=='7')]='Stem'
pbmc@meta.data$newall[which(pbmc@meta.data$clust %in% c('16','22'))]='TA'
pbmc@meta.data$newall[which(pbmc@meta.data$clust %in% c('6','15','8','24'))]='ImmuneCell'
pbmc@meta.data$newall[which(!pbmc@meta.data$newall %in% c('Stem','TA','ImmuneCell'))]='Enterocyte'

DimPlot(pbmc, group.by='newall',label=T)
VlnPlot(pbmc,group.by='newall',features=c('Olfm4'))
VlnPlot(pbmc,group.by='newall',features=c('Mki67'))




DimPlot(pbmc, group.by='newall',label=T)
#source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
source('scRef.R')

USED=which(pbmc@meta.data$newall %in% c('Enterocyte'))
exp_sc=as.matrix(pbmc@assays$RNA@data)
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
VlnPlot(pbmc,group.by='newall',features=c('Olfm4'))
VlnPlot(pbmc,group.by='newall',features=c('Mki67'))
dev.off()




OUT=t(table(pbmc@meta.data$batch,pbmc@meta.data$newall))
write.table(OUT,file='BSTAT.txt',sep='\t',row.names=T,col.names=T,quote=F)





STEM=which(pbmc@meta.data$newall %in% c('Stem','TA'))
EXP=as.matrix(pbmc@assays$RNA@data[,STEM])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(pbmc@meta.data$batch[STEM]))
OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'PT.cls',sep=' ',quote=F,row.names=F,col.names=F )


STEM=which(pbmc@meta.data$newall %in% c('Stem'))
EXP=as.matrix(pbmc@assays$RNA@data[,STEM])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(pbmc@meta.data$batch[STEM]))
OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'STEMEXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'STEMPT.cls',sep=' ',quote=F,row.names=F,col.names=F )







USED_GENE1=c('CDC20','GADD45B','SLIT2','TSPAN3','BIRC5','PHGDH','EMP2','FLNA','GLS','SH2D4A','MARCKS','GGH','SHCBP1',
             'STMN1','DUT','TK1','ETV5','TOP2A','CENPF','HMMR')



EXP=as.matrix(pbmc@assays$RNA@data[,STEM])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
META=pbmc@meta.data[STEM,]
COL=rep('blue',ncol(EXP))
COL[which(META$batch=='CDC42KO')]='red'

#EXP=log(EXP+1,10)
SEXP=EXP
#SEXP=t(apply(EXP,1,scale))
SEXP=apply(SEXP,2,pnorm)
colnames(SEXP)=colnames(EXP)
rownames(SEXP)=rownames(SEXP)

library('gplots')
heatmap.2(EXP[which(toupper(rownames(EXP)) %in% USED_GENE1),],
          ColSideColors=COL,labCol='',scale=c("none"),dendrogram='none',
          Rowv=T,Colv=F,trace='none',col=colorRampPalette(c('blue','yellow','red1','red3')),margins=c(5,5))


