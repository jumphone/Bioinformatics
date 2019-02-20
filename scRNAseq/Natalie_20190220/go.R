
library('Seurat')
source('scRef.R')

load('Seurat_EXP_cluster.Robj')

source('scRef.R')
pbmc=EXP_cluster
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

tmp0=readRDS('Inj9dBeads_Ref.RDS')
tmp0=tmp0[,which(!colnames(tmp0) %in% c('Mesenchymal.cells','Dividing.mesenchymal.cells'))]
colnames(tmp0)=paste0('Inj9d.',colnames(tmp0))

tmp1=readRDS('Inj9dBeadsMesenchymal_Ref.RDS')
tmp2=readRDS('UninjMesenchymal_Ref.RDS')

colnames(tmp1)=paste0('Inj9dMesenchymal.',colnames(tmp1))
colnames(tmp2)=paste0('UninjMesenchymal.',colnames(tmp2))

exp_ref_mat=.simple_combine(tmp0, tmp1)$combine
exp_ref_mat=.simple_combine(exp_ref_mat, tmp2)$combine

#########################
out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=6, print_step=10)
tag=.get_tag_max(out)
LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)
out=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=6, print_step=10)
tag=.get_tag_max(out)


pbmc@meta.data$injmech=tag[,2]


write.table(tag,file='injmec_tag.txt',quote=F,row.names=F,col.names=T,sep='\t') 

pdf('INJMECH.pdf',width=25,height=20)
TSNEPlot(object = pbmc, group.by ='injmech', do.label=TRUE,label.size=2)
TSNEPlot(object = pbmc, group.by ='injmech', do.label=F)
dev.off()


