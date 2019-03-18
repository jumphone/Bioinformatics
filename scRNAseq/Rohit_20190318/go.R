load('GliomaSingleCell.RData')


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

mybeer=MBEER(DATA, BATCH, MAXBATCH="", CNUM=10, PCNUM=20, CPU=2, SEED=1 )

PCUSE <- which(mybeer$cor> 0.6  & mybeer$fdr<0.05)

pbmc=mybeer$seurat

pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)  


saveRDS(pbmc,'pbmc.RDS')
saveRDS(mybeer,'mybeer.RDS')

