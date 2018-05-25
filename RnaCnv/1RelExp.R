TPM_MATRIX="NSC_merge_fpkm_named.txt.pure"

suppressPackageStartupMessages(library(edgeR))
scale_center=function(x){ return(scale(x,center=T,scale=F)) }
raw_exp_data=read.delim(TPM_MATRIX, sep='\t', header=T, row.names=1)

raw_exp_data=as.matrix(raw_exp_data)
SUM=apply(raw_exp_data, 2,sum)
i=1
while(i<=length(SUM)){
raw_exp_data[,i]=raw_exp_data[,i]/SUM[i]
i=i+1}


tmm=edgeR::calcNormFactors(raw_exp_data)
exptmm=edgeR::cpm(raw_exp_data, lib.size = tmm * colSums(raw_exp_data))

#logexp <- exptmm #log2(exptmm + 1)
logexp <- log2(exptmm + 1)

agg=apply(logexp,1,sum)
#gene_analyzed= which( agg > quantile(agg,0.0))

allvar=apply(logexp,1,var)
gene_analyzed= which( agg > quantile(agg,0.05) & agg < quantile(agg,0.95) &  allvar > quantile(allvar,0.05) & allvar < quantile(allvar,0.95))


exp_data=exptmm[gene_analyzed,]
REL_EXP=t(apply(exp_data,1,scale_center))
colnames(REL_EXP)=colnames(raw_exp_data)
write.table(REL_EXP,'REL_EXP.txt',sep='\t',row.names=T,col.names=T,quote=F)

