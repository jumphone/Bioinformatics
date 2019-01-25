source('scRef.R')
library(Seurat)
exp_ref_mat=readRDS('./Sup4_projection/exp_ref_mat.RDS')
ref_tag=readRDS('./Sup4_projection/ref_tag.RDS')
ref_vec=readRDS('./Sup4_projection/ref_vec.RDS')
TAB=read.table('Sup1_CellTypeAndClust.txt',sep='\t',row.names=1,header=T)
LocalRef=read.table('Sup3_ClustRef_human.txt',sep='\t',row.names=1,header=T)

############
NUMBER=.check_pos(LocalRef)
LocalRef=.trim_pos(LocalRef, min(NUMBER))
#############

CPU=8

########################################
########################################
TAG='R4'
################
scdata=readRDS(paste0(TAG,'_umap.RDS'))
########################################
COL=c()
i=1
while(i <=length(scdata@ident)){
    this_col=which(colnames(scdata@raw.data)==names(scdata@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(scdata@raw.data)[,COL]
########################################
scref_out=SCREF(exp_sc_mat, LocalRef,CPU=CPU)
sc_tag=scref_out$tag2
##########
scdata@meta.data$scref=sc_tag[,2]
pdf(paste0('./Sup4_projection/',TAG,'_tsne.pdf'),width=7,height=7)
TSNEPlot(scdata,group.by='scref',do.label=TRUE,pt.size=0.5)
dev.off()
################
this_tab=table(sc_tag[,2])
this_lab=c()
for(one in names(this_tab)){ 
    tmp_lab=colnames(TAB)[which(TAB[which(rownames(TAB)==one),]>0)]
    this_lab=c(this_lab, tmp_lab)
}
################
this_stat=cbind(names(this_tab),this_lab,this_tab)
colnames(this_stat)=c('cluster','cell_type','cell_num')
write.table(this_stat,paste0('./Sup4_projection/',TAG,'_stat.txt'),sep='\t',quote=F,row.names=F,col.names=F)
################
out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=50, 
        random_seed=123, CPU=CPU, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(scref_out,paste0('./Sup4_projection/',TAG,'_scref_out.RDS'))
saveRDS(out,paste0('./Sup4_projection/',TAG,'_out.RDS'))
################
pdf(paste0('./Sup4_projection/',TAG,'_proj.pdf'),width=7,height=7)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.3)
dev.off()

########################################
########################################

TAG='R18058913'
################
scdata=readRDS(paste0(TAG,'_umap.RDS'))
########################################
COL=c()
i=1
while(i <=length(scdata@ident)){
    this_col=which(colnames(scdata@raw.data)==names(scdata@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(scdata@raw.data)[,COL]
########################################
scref_out=SCREF(exp_sc_mat, LocalRef,CPU=CPU)
sc_tag=scref_out$tag2
##########
scdata@meta.data$scref=sc_tag[,2]
pdf(paste0('./Sup4_projection/',TAG,'_tsne.pdf'),width=7,height=7)
TSNEPlot(scdata,group.by='scref',do.label=TRUE,pt.size=0.5)
dev.off()
################
this_tab=table(sc_tag[,2])
this_lab=c()
for(one in names(this_tab)){ 
    tmp_lab=colnames(TAB)[which(TAB[which(rownames(TAB)==one),]>0)]
    this_lab=c(this_lab, tmp_lab)
}
################
this_stat=cbind(names(this_tab),this_lab,this_tab)
colnames(this_stat)=c('cluster','cell_type','cell_num')
write.table(this_stat,paste0('./Sup4_projection/',TAG,'_stat.txt'),sep='\t',quote=F,row.names=F,col.names=F)
################
out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=50, 
        random_seed=123, CPU=CPU, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(scref_out,paste0('./Sup4_projection/',TAG,'_scref_out.RDS'))
saveRDS(out,paste0('./Sup4_projection/',TAG,'_out.RDS'))
################
pdf(paste0('./Sup4_projection/',TAG,'_proj.pdf'),width=7,height=7)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.3)
dev.off()

########################################
########################################

TAG='R18058914'
################
scdata=readRDS(paste0(TAG,'_umap.RDS'))
########################################
COL=c()
i=1
while(i <=length(scdata@ident)){
    this_col=which(colnames(scdata@raw.data)==names(scdata@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(scdata@raw.data)[,COL]
########################################
scref_out=SCREF(exp_sc_mat, LocalRef,CPU=CPU)
sc_tag=scref_out$tag2
##########
scdata@meta.data$scref=sc_tag[,2]
pdf(paste0('./Sup4_projection/',TAG,'_tsne.pdf'),width=7,height=7)
TSNEPlot(scdata,group.by='scref',do.label=TRUE,pt.size=0.5)
dev.off()
################
this_tab=table(sc_tag[,2])
this_lab=c()
for(one in names(this_tab)){ 
    tmp_lab=colnames(TAB)[which(TAB[which(rownames(TAB)==one),]>0)]
    this_lab=c(this_lab, tmp_lab)
}
################
this_stat=cbind(names(this_tab),this_lab,this_tab)
colnames(this_stat)=c('cluster','cell_type','cell_num')
write.table(this_stat,paste0('./Sup4_projection/',TAG,'_stat.txt'),sep='\t',quote=F,row.names=F,col.names=F)
################
out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=50, 
        random_seed=123, CPU=CPU, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(scref_out,paste0('./Sup4_projection/',TAG,'_scref_out.RDS'))
saveRDS(out,paste0('./Sup4_projection/',TAG,'_out.RDS'))
################
pdf(paste0('./Sup4_projection/',TAG,'_proj.pdf'),width=7,height=7)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.3)
dev.off()

########################################
########################################

TAG='R18059655'
################
scdata=readRDS(paste0(TAG,'_umap.RDS'))
########################################
COL=c()
i=1
while(i <=length(scdata@ident)){
    this_col=which(colnames(scdata@raw.data)==names(scdata@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(scdata@raw.data)[,COL]
########################################
scref_out=SCREF(exp_sc_mat, LocalRef,CPU=CPU)
sc_tag=scref_out$tag2
##########
scdata@meta.data$scref=sc_tag[,2]
pdf(paste0('./Sup4_projection/',TAG,'_tsne.pdf'),width=7,height=7)
TSNEPlot(scdata,group.by='scref',do.label=TRUE,pt.size=0.5)
dev.off()
################
this_tab=table(sc_tag[,2])
this_lab=c()
for(one in names(this_tab)){ 
    tmp_lab=colnames(TAB)[which(TAB[which(rownames(TAB)==one),]>0)]
    this_lab=c(this_lab, tmp_lab)
}
################
this_stat=cbind(names(this_tab),this_lab,this_tab)
colnames(this_stat)=c('cluster','cell_type','cell_num')
write.table(this_stat,paste0('./Sup4_projection/',TAG,'_stat.txt'),sep='\t',quote=F,row.names=F,col.names=F)
################
out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=50, 
        random_seed=123, CPU=CPU, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(scref_out,paste0('./Sup4_projection/',TAG,'_scref_out.RDS'))
saveRDS(out,paste0('./Sup4_projection/',TAG,'_out.RDS'))
################
pdf(paste0('./Sup4_projection/',TAG,'_proj.pdf'),width=7,height=7)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.3)
dev.off()

########################################
########################################

TAG='BT310'
################
scdata=readRDS(paste0(TAG,'_umap.RDS'))
########################################
COL=c()
i=1
while(i <=length(scdata@ident)){
    this_col=which(colnames(scdata@raw.data)==names(scdata@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(scdata@raw.data)[,COL]
########################################
scref_out=SCREF(exp_sc_mat, LocalRef,CPU=CPU)
sc_tag=scref_out$tag2
##########
scdata@meta.data$scref=sc_tag[,2]
pdf(paste0('./Sup4_projection/',TAG,'_tsne.pdf'),width=7,height=7)
TSNEPlot(scdata,group.by='scref',do.label=TRUE,pt.size=0.5)
dev.off()
################
this_tab=table(sc_tag[,2])
this_lab=c()
for(one in names(this_tab)){ 
    tmp_lab=colnames(TAB)[which(TAB[which(rownames(TAB)==one),]>0)]
    this_lab=c(this_lab, tmp_lab)
}
################
this_stat=cbind(names(this_tab),this_lab,this_tab)
colnames(this_stat)=c('cluster','cell_type','cell_num')
write.table(this_stat,paste0('./Sup4_projection/',TAG,'_stat.txt'),sep='\t',quote=F,row.names=F,col.names=F)
################
out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=50, 
        random_seed=123, CPU=CPU, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(scref_out,paste0('./Sup4_projection/',TAG,'_scref_out.RDS'))
saveRDS(out,paste0('./Sup4_projection/',TAG,'_out.RDS'))
################
pdf(paste0('./Sup4_projection/',TAG,'_proj.pdf'),width=7,height=7)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.3)
dev.off()


########################################
########################################

TAG='BT311'
################
scdata=readRDS(paste0(TAG,'_umap.RDS'))
########################################
COL=c()
i=1
while(i <=length(scdata@ident)){
    this_col=which(colnames(scdata@raw.data)==names(scdata@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(scdata@raw.data)[,COL]
########################################
scref_out=SCREF(exp_sc_mat, LocalRef,CPU=CPU)
sc_tag=scref_out$tag2
##########
scdata@meta.data$scref=sc_tag[,2]
pdf(paste0('./Sup4_projection/',TAG,'_tsne.pdf'),width=7,height=7)
TSNEPlot(scdata,group.by='scref',do.label=TRUE,pt.size=0.5)
dev.off()
################
this_tab=table(sc_tag[,2])
this_lab=c()
for(one in names(this_tab)){ 
    tmp_lab=colnames(TAB)[which(TAB[which(rownames(TAB)==one),]>0)]
    this_lab=c(this_lab, tmp_lab)
}
################
this_stat=cbind(names(this_tab),this_lab,this_tab)
colnames(this_stat)=c('cluster','cell_type','cell_num')
write.table(this_stat,paste0('./Sup4_projection/',TAG,'_stat.txt'),sep='\t',quote=F,row.names=F,col.names=F)
################
out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=50, 
        random_seed=123, CPU=CPU, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(scref_out,paste0('./Sup4_projection/',TAG,'_scref_out.RDS'))
saveRDS(out,paste0('./Sup4_projection/',TAG,'_out.RDS'))
################
pdf(paste0('./Sup4_projection/',TAG,'_proj.pdf'),width=7,height=7)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.3)
dev.off()
