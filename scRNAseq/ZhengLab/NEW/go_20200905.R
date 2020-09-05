setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

library(Seurat)
pbmc=readRDS('../pbmc.final.RDS')

wt.pbmc=readRDS('new.wt.pbmc.rds')
ko.pbmc=readRDS('new.ko.pbmc.rds')

HQRef=readRDS('HQRef.rds')
ref_tag=cbind(colnames(pbmc),as.character(pbmc@meta.data$celltype))

exp_ref_mat=as.matrix(pbmc@assays$RNA@data)
exp_sc_mat=as.matrix(wt.pbmc@assays$RNA@data)

out=.get_cor(exp_sc_mat, HQRef, method='spearman',CPU=4, print_step=10)
sc_tag=.get_tag_max(out)

ref_vec=pbmc@reductions$umap@cell.embeddings
saveRDS(ref_vec, 'ref_vec.rds')



pbmc=readRDS('../pbmc.final.RDS')
VG=VariableFeatures(pbmc)
saveRDS(VG, file='VG.rds')

exp_ref_mat[which(rownames(exp_ref_mat) %in% VariableFeatures(pbmc)),]

######################################
setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

library(Seurat)
exp_ref_mat=readRDS('exp_ref_mat.rds')
ref_vec=readRDS('ref_vec.rds')
ref_tag=readRDS('ref_tag.rds')
HQRef=readRDS('HQRef.rds')
VG=readRDS('VG.rds')
wt.pbmc=readRDS('new.wt.pbmc.rds')
ko.pbmc=readRDS('new.ko.pbmc.rds')

exp_sc_mat=as.matrix(wt.pbmc@assays$RNA@data)
sc_tag=readRDS('wt.sc_tag.rds')

out =.vec_projection(exp_sc_mat[which(rownames(exp_sc_mat) %in% VG),], sc_tag, exp_ref_mat[which(rownames(exp_ref_mat) %in% VG),], ref_tag, ref_vec, 
        method='spearman', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, CPU=2, print_step=10)







