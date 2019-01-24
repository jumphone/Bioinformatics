library(Seurat)
library(dplyr)

pbmc.data <- readRDS('cerebullum_dev.RDS')
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
      
length(x = pbmc@var.genes)
#4602

pbmc <- ScaleData(object = pbmc, genes.use =pbmc@var.genes, vars.to.regress = c("nUMI", "percent.mito"),num.cores =4,do.par=TRUE)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

meta.data=read.table('cerebellum_cell_metadata.tsv',sep='\t',header=T)

pbmc@dr$tsne@cell.embeddings[,1]=as.numeric(meta.data[,1])
pbmc@dr$tsne@cell.embeddings[,2]=as.numeric(meta.data[,2])
pbmc@meta.data$tag=as.character(meta.data[,15])
pbmc@meta.data$dev=as.character(meta.data[,4])
pbmc@meta.data$clust=meta.data[,14]

saveRDS(pbmc,'cb_seurat.RDS')
########################################

pbmc=readRDS('cb_seurat.RDS')

TAB=table(pbmc@meta.data$clust,pbmc@meta.data$tag)
write.table(TAB,file='Sup1_CellTypeAndClust.txt',sep='\t',quote=F,row.names=T,col.names=T)

pdf('Sup2_CellTypeAndClust.pdf',width=20,height=18)
TSNEPlot(pbmc,group.by='tag',do.label=TRUE)
TSNEPlot(pbmc,group.by='clust',do.label=TRUE)
dev.off()

source('scRef.R')
ref_vec=pbmc@dr$tsne@cell.embeddings
ref_tag=cbind(names(pbmc@ident), as.character(pbmc@meta.data$clust))    
exp_ref_mat=as.matrix(pbmc@raw.data)
rownames(exp_ref_mat)=toupper(rownames(exp_ref_mat))

#####remove dup gene name###
tmp_tab=table(rownames(exp_ref_mat))
tmp_gene=names(which(tmp_tab>1))
exp_ref_mat=exp_ref_mat[which(!rownames(exp_ref_mat) %in% tmp_gene),]
########################################

LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )  
write.table(LocalRef,'Sup3_ClustRef.txt',sep='\t',quote=F,row.names=T,col.names=T)
saveRDS(exp_ref_mat,'./Sup4_projection/exp_ref_mat.RDS')
saveRDS(ref_tag,'./Sup4_projection/ref_tag.RDS')
################################################################################
################################################################################
########################################
source('scRef.R')
exp_ref_mat=readRDS('./Sup4_projection/exp_ref_mat.RDS')
ref_tag=readRDS('./Sup4_projection/ref_tag.RDS')
TAB=read.table('Sup1_CellTypeAndClust.txt',sep='\t',row.names=1,header=T)
LocalRef=read.table('Sup3_ClustRef.txt',sep='\t',row.names=1,header=T)
########################################
TAG='R4'
################
scdata=readRDS(paste0(TAG,'_umap.RDS'))
exp_sc_mat=as.matrix(scdata@raw.data) 
sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2
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
        method='kendall', nearest_cell=1, alpha=0.5, random_size=300, 
        random_seed=123, CPU=4, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(sc_tag,paste0('./Sup4_projection/',TAG,'_sc_tag.RDS'))
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
exp_sc_mat=as.matrix(scdata@raw.data) 
sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2
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
        method='kendall', nearest_cell=1, alpha=0.5, random_size=300, 
        random_seed=123, CPU=4, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(sc_tag,paste0('./Sup4_projection/',TAG,'_sc_tag.RDS'))
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
exp_sc_mat=as.matrix(scdata@raw.data) 
sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2
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
        method='kendall', nearest_cell=1, alpha=0.5, random_size=300, 
        random_seed=123, CPU=4, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(sc_tag,paste0('./Sup4_projection/',TAG,'_sc_tag.RDS'))
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
exp_sc_mat=as.matrix(scdata@raw.data) 
sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2
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
        method='kendall', nearest_cell=1, alpha=0.5, random_size=300, 
        random_seed=123, CPU=4, print_step=10)
################
saveRDS(exp_sc_mat,paste0('./Sup4_projection/',TAG,'_exp_sc_mat.RDS'))
saveRDS(sc_tag,paste0('./Sup4_projection/',TAG,'_sc_tag.RDS'))
saveRDS(out,paste0('./Sup4_projection/',TAG,'_out.RDS'))
################
pdf(paste0('./Sup4_projection/',TAG,'_proj.pdf'),width=7,height=7)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.3)
dev.off()







