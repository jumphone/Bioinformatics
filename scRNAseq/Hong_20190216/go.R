library(Seurat)
library(dplyr)

#smart
data_6678 <- read.table('ref_6678_pure.txt',header=T,row.names=1,sep='\t',check.names=F)
#############
gene_info=read.table('Homo_sapiens.GRCh38.87.chr.gtf.combined.pc.bed',header=F,sep='\t',check.names=F)
used=which(rownames(data_6678) %in% gene_info[,5])
data_6678=data_6678[used,]
i=1
while(i<=nrow(data_6678)){
    this_id=which(gene_info[,5]==rownames(data_6678)[i])
    data_6678[i,]=data_6678[i,]/(gene_info[this_id,3]-gene_info[this_id,2])
    i=i+1
    if(i%%10==1){print(i)}}
################
data_6678[1:3,1:3]
colnames(data_6678)=paste0('cluster',colnames(data_6678))
data_6678[1:3,1:3]

#drop
data_6701 <- read.table('ref_6701_pure.txt',header=T,row.names=1,sep='\t',check.names=F)
data_6701[1:3,1:3]
colnames(data_6701)=paste0('cluster',colnames(data_6701))
data_6701[1:3,1:3]

source('scRef.R')
n6678=.check_pos(data_6678)
n6701=.check_pos(data_6701)

data_6678_trim=.trim_pos(data_6678, min(n6678))  
data_6701_trim=.trim_pos(data_6701, min(n6701))  

exp_sc_mat=data_6678_trim
exp_ref_mat=data_6701_trim
out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
tag=.get_tag_max(out)


exp_sc_mat=data_6701 #drop
exp_ref_mat=data_6678 #smart
out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
tag=.get_tag_max(out)



