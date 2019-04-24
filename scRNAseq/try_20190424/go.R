library(Seurat)
library(dplyr)
library(Matrix)
pbmc.data <- read.table("MGH54_mat.txt", sep='\t',row.names=1, header=T)
LR=read.table('RL.txt',header=T,sep='\t')

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')


ONE=.data2one(pbmc.data, rownames(pbmc.data), CPU=4, PCNUM=50, SEED=123,  PP=30)

WINDOW=50

ORDER=order(ONE)
LENGTH=length(ONE)
BIN=c()
i=1
while(i<=LENGTH/WINDOW  ){
this_index=which((i-1)*WINDOW< ORDER & i*WINDOW>=ORDER)
BIN=cbind(BIN,this_index)
i=i+1
}

EXP=pbmc.data
GENE=rownames(EXP)



################################################

ALL=c(as.character(LR[,1]),as.character(LR[,2]) )

permu_gene_index=which(GENE %in% ALL)

set.seed(123)
TIME=1000
MEAN=EXP[permu_gene_index,c(1:TIME)]*0
colnames(MEAN)=as.character(c(1:TIME))
i=1
while(i<=TIME){
this_index=sample(ncol(EXP),WINDOW)
this_mean=apply(EXP[permu_gene_index,this_index],1,mean)
MEAN[,i]=this_mean
i=i+1
if(i%%100==1){print(i)}}

saveRDS(MEAN,file=paste0('MEAN.RDS' ))

################################################
EXP_LR=EXP[permu_gene_index,]


ECDF=c()
j=1
while(j<=nrow(MEAN)){
ECDF=c(ECDF,ecdf(as.numeric(MEAN[j,])))
j=j+1}


PMAT = MEAN[,c(1:ncol(BIN))]*0
i=1
while(i<=ncol(BIN)){
this_bin_index=BIN[,i]
  
this_bin_mean_exp=apply(EXP_LR[,this_bin_index],1,mean)  
this_p_list=c()
j=1
while(j<=length(this_bin_mean_exp)){
  this_p=-log(1+1/TIME-ECDF[[j]](this_bin_mean_exp[j]),10)
  this_p_list=c(this_p_list,this_p)
  j=j+1
}
PMAT[,i]=this_p_list
print(i)
i=i+1
}
saveRDS(PMAT,file=paste0('PMAT.RDS' ))

################################################
CMAT=PMAT[c(1:ncol(PMAT)),]*0
rownames(CMAT)=colnames(CMAT)
rownames(CMAT)=paste0('L_',rownames(CMAT))
colnames(CMAT)=paste0('R_',colnames(CMAT))

i=1
while(i<=nrow(LR)){

this_l=LR[i,1]
this_r=LR[i,2]
if(this_l %in% GENE & this_r %in% GENE){
    this_l_index=which(rownames(PMAT)==this_l)
    this_r_index=which(rownames(PMAT)==this_r)
    this_l_bin_index=1
    while(this_l_bin_index<=nrow(CMAT)){
       this_r_bin_index=1
       while(this_r_bin_index<=ncol(CMAT)){
           CMAT[this_l_bin_index,this_r_bin_index]=CMAT[this_l_bin_index,this_r_bin_index]+ 
               PMAT[this_l_index,this_l_bin_index] + PMAT[this_r_index,this_r_bin_index] 
           this_r_bin_index=this_r_bin_index+1
           }      
       this_l_bin_index=this_l_bin_index+1
       } 
    }
if(i%%10==1){print(i)}
i=i+1}

CMAT=as.matrix(CMAT)

heatmap(CMAT,scale='none')

library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))


BIN_FLAG=rep(NA,ncol(EXP))
i=1
while(i<=ncol(BIN)){
BIN_FLAG[BIN[,i]]=i
i=i+1
}

vis_gene='PDGFRA'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='SOX2'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='OLIG2'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='GFAP'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='CD83'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

plot(as.numeric(EXP[which(GENE==vis_gene),])~jitter(BIN_FLAG))

library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))










