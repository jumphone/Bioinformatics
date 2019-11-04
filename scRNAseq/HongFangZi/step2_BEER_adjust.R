setwd('/data/tianlab/zhangfeng/HFZ')
print(getwd())

source('BEER.R')
#####################
mybeer=readRDS(file='mybeer.RDS')





pbmc <- mybeer$seurat
PCUSE=mybeer$select
pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat

#umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
#pbmc@reductions$umap@cell.embeddings=umap

pdf('~/Downloads/HFZ2.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
dev.off()

saveRDS(pbmc, file='pbmc.enh.RDS')
