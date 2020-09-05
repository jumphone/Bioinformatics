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

out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='spearman', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, CPU=2, print_step=10)


