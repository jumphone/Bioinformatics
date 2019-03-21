library(Seurat)
load('Lats12SN2samples.RData')
DimPlot(GFPSN_merged, reduction.use='tsne', group.by='orig.ident')

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


pbmc=GFPSN_merged
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]



###MCA

exp_ref_mat=read.table('combined_ref.txt',sep='\t',row.names=1,header=T)

tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
pbmc@meta.data$mca=tag[,2]
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='mca')

###BrainDev

exp_ref_mat=read.table('exp_ref_mat_mouse_brain_dev.txt',header=T,row.names=1,sep='\t',check.name=F)
tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
pbmc@meta.data$dev=tag[,2]
COLOR=heat.colors(n=length(table(pbmc@meta.data$scref))+2)
TSNEPlot(object = pbmc, colors.use=COLOR, group.by ='dev')

###Inj_bulk

exp_ref_mat=read.table('injury_ref_mouse.txt',header=T,row.names=1,sep='\t',check.name=F)
tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
pbmc@meta.data$injbulk=tag[,2]
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='injbulk')

###Inj_single

tmp0=readRDS('Inj9dBeads_Ref.RDS')
tmp0=tmp0[,which(!colnames(tmp0) %in% c('Mesenchymal.cells','Dividing.mesenchymal.cells'))]
colnames(tmp0)=paste0('Inj9d.',colnames(tmp0))

tmp1=readRDS('Inj9dBeadsMesenchymal_Ref.RDS')
tmp2=readRDS('UninjMesenchymal_Ref.RDS')

colnames(tmp1)=paste0('Inj9dMesenchymal.',colnames(tmp1))
colnames(tmp2)=paste0('UninjMesenchymal.',colnames(tmp2))

exp_ref_mat=.simple_combine(tmp0, tmp1)$combine
exp_ref_mat=.simple_combine(exp_ref_mat, tmp2)$combine

tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
pbmc@meta.data$injsc=tag[,2]
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='injsc')




