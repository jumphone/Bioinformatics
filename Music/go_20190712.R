D1=readRDS(file='D1.RDS')
D2=readRDS(file='D2.RDS')

T1=read.table('E14_combined_matrix_ClusterAnnotations.txt',sep='\t',row.names=1,header=T)
T1[,1]=as.character(T1[,1])
AT1=c()
i=1
while(i<=ncol(D1)){
this_t='NA'
if(colnames(D1)[i] %in% rownames(T1)){
    this_t=T1[which(rownames(T1)==colnames(D1)[i]),1]    
    }
AT1=c(AT1,this_t)
i=i+1
}


T2=read.table('P0_combined_matrix_ClusterAnnotations.txt',sep='\t',row.names=1,header=T)
T2[,1]=as.character(T2[,1])
AT2=c()
i=1
while(i<=ncol(D2)){
this_t='NA'
if(colnames(D2)[i] %in% rownames(T2)){
    this_t=T2[which(rownames(T2)==colnames(D2)[i]),1]    
    }
AT2=c(AT2,this_t)
i=i+1
}



source('Delia.R')
COM=.simple_combine(D1,D2)$combine
COM=as.matrix(COM)
CELLTYPE=c(AT1,AT2)

REF=.generate_ref(COM, CELLTYPE)
REF=log(REF+1,10)


EXP=read.table('gene_rep_matix_anno_diff.txt.pure',sep='\t',row.names=1,header=T)
EXP=as.matrix(EXP)
EXP=log(EXP+1,10)




mydelia <- Delia(EXP, REF, COMBAT=TRUE, METHOD='rlm')  




show_ratio_coef <- mydelia$coef
#show_ratio_coef = t(apply(show_ratio_coef,1,scale))
rownames(show_ratio_coef )=rownames(mydelia$coef)
colnames(show_ratio_coef )=colnames(mydelia$coef)

library('gplots')
 


as.numeric(show_ratio_coef)

DATA=mydelia$coef
OUT=matrix(scale(as.numeric(DATA)), ncol=ncol(DATA),nrow=nrow(DATA))
#OUT=pnorm(OUT)
rownames(OUT)=rownames(DATA)
colnames(OUT)=colnames(DATA)



show_ratio_coef=OUT
heatmap.2(t(show_ratio_coef),scale=c("none"), dendrogram='column',
    Rowv=F,Colv=T,cellnote=round(t(show_ratio_coef),2), notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
    margins=c(10,10))




show_ratio_coef <- mydelia$out

library('gplots')
 
heatmap.2(t(show_ratio_coef),scale=c("none"), dendrogram='none',
    Rowv=F,Colv=F,cellnote=round(t(show_ratio_coef),2), notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
    margins=c(10,10))



