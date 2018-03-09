RE_MATRIX="REL_EXP.sorted.txt"

WINDOW=100
scale_center=function(x){return(scale(x,center=T,scale=F))}

raw_data=read.delim(RE_MATRIX)
relative_gene_exp=as.matrix(raw_data[,c(5:length(raw_data[1,]))])
relative_gene_exp[which(relative_gene_exp > 3)]= 3
relative_gene_exp[which(relative_gene_exp < -3)]= -3
cell_analyzed_num=length(relative_gene_exp[1,])
gene_analyzed_num=length(relative_gene_exp[,1])
info_relative_gene_exp= raw_data[,c(1:4)]
CHROM=info_relative_gene_exp[,1]





ALLCNV=c()
j=1
while(j<=cell_analyzed_num){
CNV=c()
TPM=c()
begin=WINDOW/2+1
end=gene_analyzed_num-WINDOW/2
i=begin
while(i <=end){
sss=(i-WINDOW/2)
eee=(i+WINDOW/2)
if(CHROM[sss]!=CHROM[i]){sss= i}
if(CHROM[eee]!=CHROM[i]){eee= i}
CNV=c(CNV,sum(relative_gene_exp[c(sss:eee),j])/(eee-sss+1))
i=i+1}
ALLCNV=cbind(ALLCNV,CNV)
j=j+1
print(j)
}
rownames(ALLCNV) = rownames(relative_gene_exp)[begin:end]
colnames(ALLCNV) = colnames(relative_gene_exp)
info_ALLCNV = info_relative_gene_exp[begin:end,]
center_ALLCNV=apply(ALLCNV,2,scale_center)
rownames(center_ALLCNV)=rownames(ALLCNV)
colnames(center_ALLCNV)=colnames(ALLCNV)
write.table(cbind(info_ALLCNV,center_ALLCNV), file='center_ALLCNV.txt',row.names=F,col.names=T,quote=F,sep='\t')




