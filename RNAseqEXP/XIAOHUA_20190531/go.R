
#install.packages("tidyverse")
#install.packages("readxl")
library(readxl)

DATA=read_excel("REF.xls")

MDATA=as.matrix(DATA)
RN=MDATA[,1]
MDATA=MDATA[,c(2:ncol(MDATA))]
MDATA=apply(MDATA,2,as.numeric)
rownames(MDATA)=RN
MDATA[which(is.na(MDATA))]=0
rownames(MDATA)=toupper(rownames(MDATA))
MDATA[1:3,1:3]
REF=MDATA



DATA=read_excel("NEW.xlsx")
MDATA=as.matrix(DATA)
RN=MDATA[,2]
MDATA=MDATA[,c(3:ncol(MDATA))]
MDATA=apply(MDATA,2,as.numeric)
rownames(MDATA)=RN
MDATA[which(is.na(MDATA))]=0
rownames(MDATA)=toupper(rownames(MDATA))
MDATA[1:3,1:3]
NEW=MDATA

NEW=NEW[which(rownames(NEW) %in% names(which(table(rownames(NEW))==1))),]
REF=REF[which(rownames(REF) %in% names(which(table(rownames(REF))==1))),]


SUM=apply(REF,2,sum)
plot(SUM)

expmat=REF

getNorm=function(x){

   return(x/sum(x)*1000000)
    }
expmat=apply(expmat,2,getNorm)
rownames(expmat)=rownames(REF)

VAR=apply(expmat,1,var)
expmat=expmat[which(VAR>1),]


getBatch=function(x){
   return(substr(x, 1, 3))
   }

batch=apply(matrix(colnames(REF),ncol=1),1,getBatch)

library(sva)
library(limma)
pheno = data.frame(batch=as.matrix(batch))
edata = expmat
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

SUM=apply(combat_edata,2,sum)
plot(SUM)





source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

out=.get_cor(NEW, combat_edata, method='spearman',CPU=1, print_step=10)

library('gplots')
heatmap.2(out,scale=c("none"),dendrogram='both',trace='none',col=colorRampPalette(c('blue','white','red')),margins=c(5,5))


tag=.get_tag_max(out)
write.table(tag,file='LABEL.txt',quote=F,row.names=F,col.names=T,sep='\t')
