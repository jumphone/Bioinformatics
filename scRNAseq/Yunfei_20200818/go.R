
setwd('E:/Project/YunFEI_20200818/DIPG correlation')

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


REF_MAT=read.table('TCGA_unified_CORE_ClaNC840.txt',header=T,row.names=1,sep='\t')

REF_MAT=REF_MAT[,2:ncol(REF_MAT)]
TYPE=REF_MAT[1,]
REF_MAT=REF_MAT[2:nrow(REF_MAT),]
RNAME=rownames(REF_MAT)
REF_MAT=apply(REF_MAT,2,as.numeric)
rownames(REF_MAT)=RNAME

TYPE=t(TYPE)
TYPE=cbind(colnames(REF_MAT),TYPE)
LocalRef=.generate_ref(REF_MAT, TYPE, M='mean', min_cell=10)  


SC_MAT_1=read.table('DIPG_KOOLIG2-subtype.txt',header=T,row.names=1,sep='\t')



out_1=.get_cor(SC_MAT_1, LocalRef, method='spearman',CPU=4, print_step=10)


