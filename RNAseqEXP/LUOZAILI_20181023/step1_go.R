a=read.table('2018-2-24 600 samples MB genes_matrix_symbol.txt',sep='\t')

b=read.table('2018-2-24 600 samples MB genes_matrix_symbol.txt',sep='\t',skip=3,row.names=1)

ori_groups = as.character(t(a[2,2:length(a[1,])]))

table(ori_groups)

G3=which(ori_groups %in% c('Group3_alpha','Group3_beta','Group3_gamma') )
G4=which(ori_groups %in% c('Group4_alpha','Group4_beta','Group4_gamma') )
SHH=which(ori_groups %in% c('SHH_alpha','SHH_beta','SHH_delta','SHH_gamma') )
WNT=which(ori_groups %in% c('WNT_alpha','WNT_beta','WNT_gamma') )

length(G3)+length(G4)+length(SHH)+length(WNT)
length(ori_groups)


G3mean =  as.matrix(apply(b[,G3],1,mean))
G4mean =  as.matrix(apply(b[,G4],1,mean))
SHHmean =  as.matrix(apply(b[,SHH],1,mean))
WNTmean =  as.matrix(apply(b[,WNT],1,mean))

ALL=cbind(G3mean,G4mean,SHHmean,WNTmean)
colnames(ALL)=c('G3','G4','SHH','WNT')

c=read.table('NSC-KOCTD.txt.pure',row.names=1,header=T,check.names=F)

exp_sc_mat=c
exp_ref_mat=ALL
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)    
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
colname_sc=colnames(exp_sc_mat)
colname_ref=colnames(exp_ref_mat)

exp_all=cbind(exp_sc_mat,exp_ref_mat)

exp_all[1:10,1:10]

pca=princomp(exp_all)
PC1=pca$loadings[,1]
PC2=pca$loadings[,2]

LABEL=rownames(pca$loadings)
pdf('PCA.pdf',width=7,height=7)
plot(PC1,PC2,pch=16, cex=1)
text(PC1, PC2, labels=LABEL, cex= 0.7, pos=c(1,2,3,4))
dev.off()

ALLPCC=c()
i=1
while(i <= length(colname_sc)){
	pccs=c()
	j=1
	while(j <= length(colname_ref) ){
        pcc=cor(exp_sc_mat[,i],exp_ref_mat[,j])
        pccs=c(pccs,pcc)
        j=j+1
       }
    ALLPCC=cbind(ALLPCC,pccs)
	i=i+1
}

rownames(ALLPCC)=colname_ref
colnames(ALLPCC)=colname_sc
library('gplots')
pdf('HEAT.pdf',width=7,height=7)
heatmap.2(as.matrix(ALLPCC),scale=c("none"),dendrogram='both',trace='none',col=colorRampPalette(c('blue','grey95','red')),margins=c(10,10))
dev.off()


write.table(ALLPCC,file='ALLPCC.txt',row.names=T,col.names=T,quote=F,sep='\t')


