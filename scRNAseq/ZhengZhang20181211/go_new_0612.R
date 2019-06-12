#Seurat 3.0

library(dplyr)
library(Seurat)

mybeer=readRDS('mybeer.RDS')
PCUSE <- which(mybeer$cor>quantile(mybeer$cor,0.3) )
npbmc <- mybeer$seurat
npbmc <- RunUMAP(object = npbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
DimPlot(npbmc, reduction = "umap")


pbmc=readRDS('ALL.RDS')


STEM=c("Lgr5", "Ascl2", "Slc12a2", "Axin2", "Olfm4", "Axin2")
TA=c("Mki67", "Cdk4", "Mcm5", "Mcm6", "Pcna")
QSTEM=c("Bmi1", "Hopx", "Lrig1", "Tert")



OLDTA=which(pbmc@meta.data$all %in% c('TA.Early','TA.G1','TA.G2'))

pbmc@reductions$pca@cell.embeddings[,PCUSE] 



stem_exp=apply(pbmc@assays$RNA@data[which(rownames(pbmc) %in% STEM),],2,min)
ta_exp=apply(pbmc@assays$RNA@data[which(rownames(pbmc) %in% TA),],2,min)
#qstem_exp=apply(pbmc@assays$RNA@data[which(rownames(pbmc) %in% QSTEM),],2,min)

C_STEM=which(stem_exp>0 & pbmc@meta.data$all %in% c('TA.Early','TA.G1','TA.G2') )
C_TA=which(ta_exp>0 & pbmc@meta.data$all %in% c('TA.Early','TA.G1','TA.G2') )


pbmc@meta.data$all[OLDTA]='NA'
#pbmc@meta.data$all[C_STEM]='Stem'
#pbmc@meta.data$all[C_TA]='TA'


MAT=cbind(pbmc@assays$RNA@counts[,C_STEM],pbmc@assays$RNA@counts[,C_TA])
TAG=c(rep('Stem',length(C_STEM)),rep('TA',length(C_TA)))
source('scRef.R')
REF=.generate_ref( MAT, cbind(TAG,TAG),min_cell=1)

INPUT=as.matrix(pbmc@assays$RNA@counts[,OLDTA])
out=.get_log_p_sc_given_ref(INPUT, REF, CPU=4, print_step=10)
tag=.get_tag_max(out)
pbmc@meta.data$all[OLDTA]=tag[,2]



NEWMETA=pbmc@meta.data
saveRDS(NEWMETA,file='NEWMETA.RDS')




pdf("ALL.pdf",width=10,height=5)
DimPlot(pbmc, group.by='all',label=T)
dev.off()

CDC42HET=which(pbmc@active.ident=='CDC42HET')
#CDC42HET=rep(CDC42HET,pbmc@meta.data$CDC42HET[CDC42HET])


CDC42KO=which(pbmc@active.ident=='CDC42KO')
#CDC42KO=rep(CDC42KO,pbmc@meta.data$CDC42KO[CDC42KO])



write.table(sort(table(pbmc@meta.data$all[CDC42HET])),file='CDC42HET.txt',sep='\t',quote=F,col.names=F,row.names=F)
write.table(sort(table(pbmc@meta.data$all[CDC42KO])),file='CDC42KO.txt',sep='\t',quote=F,col.names=F,row.names=F)



EXP=as.matrix(pbmc@assays$RNA@data[,OLDTA])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]

PT=t(as.character(pbmc@meta.data$batch[OLDTA]))



OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'PT.cls',sep=' ',quote=F,row.names=F,col.names=F )
