d1=read.table('4M_SN-GFP_out_gene_exon_tagged.dge.txt',sep='\t',row.names=1,header=T)
d2=read.table('p100_out_gene_exon_tagged.dge.txt',sep='\t',row.names=1,header=T)

colnames(d1)=paste0('4M_',colnames(d1))
colnames(d2)=paste0('p100_',colnames(d2))

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

mybeer <- BEER(d1, d2, CNUM=10, PCNUM=50, CPU=2, REGBATCH=TRUE)






par(mfrow=c(1,2))
plot(mybeer$cor, xlab='PCs', ylab="COR", pch=16)
plot(-log(mybeer$fdr,10), xlab='PCs', ylab='-log10(FDR)', pch=16)



PCUSE <- which(mybeer$cor> 0.5 & mybeer$fdr<0.05 & c(1:50)<=30 )

pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
#DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)



#pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast=T, check_duplicates=FALSE)
#DimPlot(pbmc, reduction.use='tsne', group.by='batch', pt.size=0.1)
#DimPlot(pbmc, reduction.use='tsne', group.by='map', pt.size=0.1)


##################

VEC=pbmc@dr$umap@cell.embeddings

# Here, we use the "dbscan" function to do clustering.
library("fpc")
set.seed(123)
df=VEC
db <- fpc::dbscan(df, eps = 0.25, MinPts = 5)
DC=db$cluster
pbmc@meta.data$DC=DC
DimPlot(pbmc, reduction.use='umap', group.by='DC', pt.size=0.5)

##############

library(dplyr)
tmp=pbmc@ident
pbmc@ident=as.factor(DC)
names(pbmc@ident)=names(tmp)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
############


save.image(file='combine_beer.RData')


pdf('BEER.pdf',width=15,height=12)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
DimPlot(pbmc, reduction.use='umap', group.by='DC', pt.size=0.5)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

########################

exp_sc_mat=d2
load('Seurat_EXP_cluster.Robj')

pbmc=EXP_cluster
ref_vec=pbmc@dr$tsne@cell.embeddings

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]
ref_tag=cbind(names(pbmc@ident),as.character(pbmc@ident))    
HQRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )   

sc_tag=SCREF(exp_sc_mat, HQRef)$tag2

out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, CPU=4, print_step=10)

########
pdf('OVERLAY.pdf',width=15,height=12)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70')
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red')
dev.off()

saveRDS(out,'overlay.RDS')
saveRDS(sc_tag,'sc_tag.RDS')


