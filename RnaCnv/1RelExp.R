TPM_MATRIX="RANDOM200TPM.txt"



suppressPackageStartupMessages(library(edgeR))
scale_center=function(x){ return(scale(x,center=T,scale=F)) }
raw_exp_data=read.delim(TPM_MATRIX, sep='\t', header=T, row.names=1)
tmm=edgeR::calcNormFactors(raw_exp_data)
exptmm=edgeR::cpm(raw_exp_data, lib.size = tmm * colSums(raw_exp_data))
#logexp <- log2(exptmm + 1)
agg=apply(exptmm,1,sum)
gene_analyzed= which( agg > quantil(agg,0.8))
exp_data=exptmm[gene_analyzed,]
REL_EXP=t(apply(exp_data,1,scale_center))
colnames(REL_EXP)=colnames(raw_exp_data)
write.table(REL_EXP,'REL_EXP.txt',sep='\t',row.names=T,col.names=T,quote=F)
