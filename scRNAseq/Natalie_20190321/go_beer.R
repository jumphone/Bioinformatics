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



pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast=T, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='tsne', group.by='batch', pt.size=0.1)
#DimPlot(pbmc, reduction.use='tsne', group.by='map', pt.size=0.1)


VEC=pbmc@dr$umap@cell.embeddings

# Here, we use the "dbscan" function to do clustering.
library("fpc")
set.seed(123)
df=VEC
db <- fpc::dbscan(df, eps = 0.5, MinPts = 5)
DC=db$cluster
pbmc@meta.data$DC=DC
DimPlot(pbmc, reduction.use='umap', group.by='DC', pt.size=0.5)

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


save(file='combine_beer.RData')
