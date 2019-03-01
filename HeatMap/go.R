source('scRef.R')
exp_ref_mat=read.table('gbm_subgroup_exp.txt',header=T,row.names=1,sep='\t')
ref_tag=read.table('gbm_subgroup_type.txt',sep='\t')

ref_tag=cbind(as.character(ref_tag[,1]),as.character(ref_tag[,1]))

exp_ref_mat=.generate_ref(exp_ref_mat,ref_tag,min_cell=1)
exp_sc_mat=read.table('DIPG_merge_fpkm.txt.pure',header=T,row.names=1,sep='\t')

out=.get_cor(exp_sc_mat,exp_ref_mat, method='spearman',CPU=4, print_step=10)
out=out[,c(1:11)]

library('gplots')

pdf('COR.pdf',width=7,height=5)
heatmap.2(out,scale=c("column"),dendrogram='none',Colv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
dev.off()

