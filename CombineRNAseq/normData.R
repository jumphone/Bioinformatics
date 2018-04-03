TPM_MATRIX='COMBINED.txt'
suppressPackageStartupMessages(library(edgeR))
raw_exp_data=read.delim(TPM_MATRIX, sep='\t', header=T, row.names=1)
tmm=edgeR::calcNormFactors(raw_exp_data)
exptmm=edgeR::cpm(raw_exp_data, lib.size = tmm * colSums(raw_exp_data))
rownames(exptmm)=rownames(raw_exp_data)
colnames(exptmm)=colnames(raw_exp_data)
write.table(exptmm,file='COMBINED_NORMED.txt',sep='\t',col.names=T,row.names=T,quote=F)


exp_data=read.table('COMBINED_NORMED.txt',header=T,row.names=1)
GENE_LIMIT=1000
GENE_NUM=c()
total_cell_number=length(exp_data[1,])
i=1
while(i <=total_cell_number){
GENE_NUM=c(GENE_NUM,length(which( exp_data[,i] > 0)))
i=i+1
print(i)}
analyzed_cell=which(GENE_NUM>=GENE_LIMIT)
exp_data_fine=exp_data[,analyzed_cell]
colnames(exp_data_fine)=colnames(exp_data)[analyzed_cell]
rownames(exp_data_fine)=rownames(exp_data)
exp_data_fine_rank=apply(-exp_data_fine, 2, rank,  ties.method='average')
colnames(exp_data_fine_rank)=colnames(exp_data_fine)
rownames(exp_data_fine_rank)=rownames(exp_data_fine)
total_gene_number=length(exp_data[,1])
average_rank = sum((GENE_LIMIT+1):total_gene_number)/(total_gene_number-GENE_LIMIT)
exp_data_fine_rank[which(exp_data_fine_rank <=GENE_LIMIT)]= GENE_LIMIT+1-exp_data_fine_rank[which(exp_data_fine_rank <=GENE_LIMIT)]
exp_data_fine_rank[which(exp_data_fine_rank >GENE_LIMIT)]=0
write.table(exp_data_fine_rank,'COMBINED_NORMED.RANK1000.txt',quote=F,sep='\t',row.names=T,col.names=T)
