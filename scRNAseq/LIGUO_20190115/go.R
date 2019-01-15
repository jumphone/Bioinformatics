source('scRef.R')

T1=read.table('10_P19_cKO1.txt.uniq.txt',sep='\t',row.names=1,header=T)
T2=read.table('11_P19_cKO2.txt.uniq.txt',sep='\t',row.names=1,header=T)
T=.simple_combine(cbind(T1,T1),cbind(T2,T2))$combine
T=T[,c(1,3)]
N=read.table('Reference_expression.txt',sep='\t',row.names=1,header=TRUE)
REF=.simple_combine(T,N)$combine
REF=REF[,c(1,3,4,5,6,7,8,9)]
library(Seurat)

T=TRUE
F=FALSE

#######################################################

load('10xEXP.Robj')
   
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
exp_ref_mat=REF 
out=.get_cor(exp_sc_mat, exp_ref_mat, method='pearson',CPU=4, print_step=10)
tag=.get_tag_max(out)
pbmc@meta.data$tag_p=tag[,2]
TOUT_P=table(pbmc@meta.data$tag_p, pbmc@ident)
write.table(TOUT_P,file='10X_pearson.txt', quote=F,row.names=TRUE,col.names=TRUE,sep='\t')

out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
tag=.get_tag_max(out)
pbmc@meta.data$tag_s=tag[,2]
TOUT_S=table(pbmc@meta.data$tag_s, pbmc@ident)
write.table(TOUT_S,file='10X_spearman.txt', quote=F,row.names=TRUE,col.names=TRUE,sep='\t')


#pdf('10X.pdf',width=7,height=5)
#TSNEPlot(pbmc,group.by='tag')
#dev.off()
scref_out=SCREF(exp_sc_mat, exp_ref_mat, min_cell=1,CPU=4, print_step=10)
scref_tag=scref_out$tag2
pbmc@meta.data$tag_scref=scref_tag[,2]
table(pbmc@meta.data$tag_scref, pbmc@ident)
TSNEPlot(pbmc,group.by='tag_scref')

#######################################################
load('dropEXP.Robj')
pbmc=EXP
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
exp_ref_mat=REF 


out=.get_cor(exp_sc_mat, exp_ref_mat, method='pearson',CPU=4, print_step=10)
tag=.get_tag_max(out)
pbmc@meta.data$tag_p=tag[,2]
TOUT_P=table(pbmc@meta.data$tag_p, pbmc@ident)
write.table(TOUT_P,file='Dropseq_pearson.txt', quote=F,row.names=TRUE,col.names=TRUE,sep='\t')

out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
tag=.get_tag_max(out)
pbmc@meta.data$tag_s=tag[,2]
TOUT_S=table(pbmc@meta.data$tag_s, pbmc@ident)
write.table(TOUT_S,file='Dropseq_spearman.txt', quote=F,row.names=TRUE,col.names=TRUE,sep='\t')



scref_out=SCREF(exp_sc_mat, exp_ref_mat, min_cell=1,CPU=4, print_step=10)
scref_tag=scref_out$tag2
pbmc@meta.data$tag_scref=scref_tag[,2]
table(pbmc@meta.data$tag_scref, pbmc@ident)
TSNEPlot(pbmc,group.by='tag_scref')
############################################################


