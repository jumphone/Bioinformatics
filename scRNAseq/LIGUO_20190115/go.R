source('scRef.R')

T1=read.table('10_P19_cKO1.txt.uniq.txt',sep='\t',row.names=1,header=T)
T2=read.table('11_P19_cKO2.txt.uniq.txt',sep='\t',row.names=1,header=T)
T=.simple_combine(cbind(T1,T1),cbind(T2,T2))$combine
T=T[,c(1,3)]
N=read.table('Reference_expression.txt',sep='\t',row.names=1,header=TRUE)
REF=.simple_combine(T,N)$combine

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
out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
tag=.get_tag_max(out)
write.table(tag,file='10X_spearman.txt',quote=F,row.names=F,col.names=T,sep='\t')
pbmc@meta.data$tag=tag[,2]
pdf('10X.pdf',width=7,height=5)
TSNEPlot(pbmc,group.by='tag')
dev.off()

#######################################################
load('DropseqEXP.Robj')
  
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
exp_ref_mat=REF 
out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
tag=.get_tag_max(out)
write.table(tag,file='Dropseq_spearman.txt',quote=F,row.names=F,col.names=T,sep='\t')
pbmc@meta.data$tag=tag[,2]
pdf('Dropseq.pdf',width=7,height=5)
TSNEPlot(pbmc,group.by='tag')
dev.off()

############################################################


