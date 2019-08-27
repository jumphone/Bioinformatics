


setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS('./pbmc1.RDS')
pbmc_zhengzhang=readRDS('pbmc_zhengzhang.RDS')

pbmc@meta.data=readRDS('pbmc1_meta.RDS')
pbmc_zhengzhang@meta.data=readRDS('pbmc_zhengzhang_meta.RDS')




DimPlot(pbmc_zhengzhang, reduction.use='umap', group.by='batch', pt.size=0.1,label=T)
DimPlot(pbmc_zhengzhang, reduction.use='umap', group.by='level1', pt.size=0.1,label=T)


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

pbmc=pbmc_zhengzhang
############################
USE.CELL=which(pbmc@meta.data$level1=='Stem')
TAG=pbmc@meta.data$batch
EXP=as.matrix(pbmc@assays$RNA@data[,USE.CELL])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(TAG[USE.CELL]))

#EXP.combat=.combat(EXP,PT[1,])
#EXP=EXP.combat
#EXP[which(EXP<0)]=0

OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'./GSEA_zhengzhang/STEM_EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'./GSEA_zhengzhang/STEM_PT.cls',sep=' ',quote=F,row.names=F,col.names=F )
##################################

############################
USE.CELL=which(pbmc@meta.data$level1=='TA')
TAG=pbmc@meta.data$batch
EXP=as.matrix(pbmc@assays$RNA@data[,USE.CELL])
VAR=apply(EXP,1,var)
EXP=EXP[which(VAR>0),]
PT=t(as.character(TAG[USE.CELL]))


#EXP.combat=.combat(EXP,PT[1,])
#EXP=EXP.combat
#EXP[which(EXP<0)]=0

OUT=cbind(toupper(rownames(EXP)),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'./GSEA_zhengzhang/TA_EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'./GSEA_zhengzhang/TA_PT.cls',sep=' ',quote=F,row.names=F,col.names=F )
##################################

