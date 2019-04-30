library(Seurat)
pbmc=readRDS('cb_seurat.RDS')


ref_ident=readRDS('pbmc_ident.rds')

tmp=as.matrix(pbmc@raw.data)
ref_exp=tmp[,c(1:length(ref_ident))]*0
colnames(ref_exp)=names(ref_ident)
i=1
while(i<=length(ref_ident)){
this_index=which(colnames(tmp)==names(ref_ident)[i])
ref_exp[,i]=tmp[,this_index]
i=i+1
if(i%%100==1){print(i)}
}


ref_tag=cbind(names(ref_ident),as.character(ref_ident))
table(ref_tag[,2])
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
LocalRef=.generate_ref(ref_exp, ref_tag, min_cell=10)  


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

rm(pbmc)
rm(ref_pbmc)
rm(ref_exp)
gc()


VAR=apply(exp_sc_mat,1,var)
used=which(VAR>=quantile(VAR,0.8))
exp_sc_mat=exp_sc_mat[used,]
gc()

out=.get_cor(exp_sc_mat, LocalRef, method='spearman',CPU=2, print_step=10)
tag=.get_tag_max(out)
saveRDS(tag,'37cluster.RDS')


rm(exp_sc_mat)
gc()
pbmc=readRDS('cb_seurat.RDS')



