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

this_l=as.character(LR[1,1])
this_r=as.character(LR[1,2])

this_l_index=which(GENE==this_l)
this_r_index=which(GENE==this_r)


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

saveRDS(MEAN,file=paste0('permutation.RDS' ))

################################################


