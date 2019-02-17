library(Seurat)
library(dplyr)


data_6678 <- read.table('ref_6678_pure.txt',header=T,row.names=1,sep='\t',check.names=F)
data_6678[1:3,1:3]
colnames(data_6678)=paste0('cluster',colnames(data_6678))
data_6678[1:3,1:3]


data_6701 <- read.table('ref_6701_pure.txt',header=T,row.names=1,sep='\t',check.names=F)
data_6701[1:3,1:3]
colnames(data_6701)=paste0('cluster',colnames(data_6701))
data_6701[1:3,1:3]

source('scRef.R')

out=.get_cor(data_6678, data_6701, method='kendall',CPU=4, print_step=10)










