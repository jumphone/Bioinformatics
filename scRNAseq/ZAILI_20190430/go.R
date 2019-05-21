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
used=which(VAR>=quantile(VAR,0.9))
exp_sc_mat=exp_sc_mat[used,]
gc()

out=.get_cor(exp_sc_mat, LocalRef, method='spearman',CPU=1, print_step=10)
tag=.get_tag_max(out)
saveRDS(tag,'37cluster.RDS')


rm(exp_sc_mat)
gc()
pbmc=readRDS('cb_seurat.RDS')
pbmc@meta.data$newtag=as.numeric(tag[,2])

tiff('PNG.tiff',width=1200,height=1000)
DimPlot(pbmc,group.by='newtag',reduction.use='tsne',do.label=T)
dev.off()

write.table(pbmc@meta.data,file='META.tsv',sep='\t',quote=F,row.names=T,col.names=T)





TAB1=table(pbmc@meta.data$dev,pbmc@meta.data$newtag)
TAB2=table(pbmc@meta.data$tag,pbmc@meta.data$newtag)


write.table(TAB1,file='TAB1.tsv',sep='\t',quote=F,row.names=T,col.names=T)
write.table(TAB2,file='TAB2.tsv',sep='\t',quote=F,row.names=T,col.names=T)


###################
library(Seurat)
pbmc=readRDS('cb_seurat.RDS')
FeaturePlot(object = pbmc, features.plot = c('Prrx2'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = pbmc, features.plot = c('Prrx2','Rgs5'), cols.use = c("grey", "blue"), reduction.use = "tsne")


########
#0521
library(Seurat)
pbmc=readRDS('cb_seurat.RDS')

ref_tag=readRDS('37cluster.RDS')
ref_map=read.table('SUB.txt',sep='\t',header=F)
i=1
while(i<=nrow(ref_map)){
this_index=which(ref_tag[,2]==as.character(ref_map[i,1]))
ref_tag[this_index,2]=as.character(ref_map[i,2])
i=i+1}


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

LocalRef=.generate_ref(ref_exp, ref_tag, min_cell=10)  

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
LocalRef=.generate_ref(exp_sc_mat, ref_tag, min_cell=10)  


saveRDS(LocalRef,file='37ref.RDS')



