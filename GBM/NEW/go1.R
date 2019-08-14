

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('/Volumes/Feng/MATT/')




#############################

library(Seurat)
library(dplyr)

pbmc=readRDS(file='newpbmc.RDS')
VF=VariableFeatures(object = pbmc)
PC1_PEAK_LOADING =pbmc@reductions$pca@feature.loadings[,1]

PC1=pbmc@reductions$pca@cell.embeddings[,1]
TAG=pbmc@meta.data$tag




D1=read.delim('matt_combined.txt',header=T,row.names=1,sep='\t')
D2=read.delim('matt_combined2.txt',header=T,row.names=1,sep='\t')


pbmc1 <- CreateSeuratObject(counts =D1, project = "D1", min.cells = 0, min.features = 0)
pbmc1 <- NormalizeData(object = pbmc1, normalization.method = "LogNormalize",  scale.factor = 10000)

pbmc2 <- CreateSeuratObject(counts =D2, project = "D2", min.cells = 0, min.features = 0)
pbmc2 <- NormalizeData(object = pbmc2, normalization.method = "LogNormalize",  scale.factor = 10000)


ND=as.matrix(pbmc@assays$RNA@data[which(rownames(pbmc) %in% VF),])
ND1=as.matrix(pbmc1@assays$RNA@data[which(rownames(pbmc1) %in% VF),])
ND2=as.matrix(pbmc2@assays$RNA@data[which(rownames(pbmc2) %in% VF),])


##########
##########
##########

COM=.simple_combine(ND,ND1)$combine
COM=.simple_combine(COM,ND2)$combine
BATCH=c(rep('ND',ncol(ND)),rep('ND1',ncol(ND1)),rep('ND2',ncol(ND2)))
COMBAT=.combat(COM,BATCH)



CND=COMBAT[,which(BATCH=='ND')]
CND1=COMBAT[,which(BATCH=='ND1')]
CND2=COMBAT[,which(BATCH=='ND2')]



source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')



#######################


exp_sc_mat=CND1
exp_ref_mat=CND
out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
tag=.get_tag_max(out)


THIS.PC1=c()
i=1
while(i<=nrow(tag)){
    THIS.PC1=c(THIS.PC1, PC1[which(names(PC1)==tag[i,2])])
    
    i=i+1}


names(THIS.PC1)=colnames(exp_sc_mat)

D1.PC1=THIS.PC1
#######################



#######################


exp_sc_mat=CND2
exp_ref_mat=CND
out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
tag=.get_tag_max(out)


THIS.PC1=c()
i=1
while(i<=nrow(tag)){
    THIS.PC1=c(THIS.PC1, PC1[which(names(PC1)==tag[i,2])])
    
    i=i+1}


names(THIS.PC1)=colnames(exp_sc_mat)

D2.PC1=THIS.PC1

#######################




SURTAG=read.table('SURTAG.txt',sep='\t',header=F)
used=which(!is.na(SURTAG[,1]))

library(survival)
library(survminer)
a=read.table('SUR.txt',header=T,sep='\t')
used=used
score=PC1[used]
SUR=SURTAG[used,2]
SURE=SURTAG[used,1]

TYPE=rep('MED',length(used))
TYPE[which(score> quantile(score,0.5) )]='High'
TYPE[which(score<= quantile(score,0.5) )]='Low'
surtime=SUR[which(TYPE!='MED')]
surevent=SURE[(which(TYPE!='MED'))]
surtype=TYPE[which(TYPE!='MED')]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
######
pdf('./FIG/1SUR.pdf',width=7,height=7)
ggsurvplot(fit, pval = TRUE)
dev.off()
######
surv_pvalue(fit)



pdf('./FIG/2PC1.pdf',width=7,height=7)
plot(sort(PC1), ylab='PC1',xlab='Rank', cex=1,pch=16)
abline(h=quantile(score,0.5),col='red',lwd=1.5)
dev.off()
######

pdf('./FIG/3MAP.pdf',width=12,height=7)
par(mar=c(15,5,5,5))
barplot(sort(c(D1.PC1,D2.PC1)),las=2,ylab="Mapped point's PC1")
dev.off()
#barplot(sort(D2.PC1),las=2,ylab='Mapped PC1')


