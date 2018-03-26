TPM_MATRIX='COMBINED.txt'
suppressPackageStartupMessages(library(edgeR))
scale_center=function(x){ return(scale(x,center=T,scale=F)) }
raw_exp_data=read.delim(TPM_MATRIX, sep='\t', header=T, row.names=1)
tmm=edgeR::calcNormFactors(raw_exp_data)
exptmm=edgeR::cpm(raw_exp_data, lib.size = tmm * colSums(raw_exp_data))

rownames(exptmm)=rownames(raw_exp_data)
colnames(exptmm)=colnames(raw_exp_data)
write.table(exptmm,file='COMBINED_NORMED.txt',sep='\t',col.names=T,row.names=T,quote=F)

