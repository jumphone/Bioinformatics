


setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS('./pbmc1.RDS')


LABEL=readRDS(file='LABEL.RDS')
pbmc@meta.data$celltype=rep(NA,ncol(pbmc))
pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='NATURE')]=LABEL
#DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=F)
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)




USED_CELL=which(pbmc@meta.data$batch %in% c('CDC42KO','CDC42HET'))
pbmc_zhengzhang <- CreateSeuratObject(counts = pbmc@assays$RNA@counts[,USED_CELL], project = "ZhengZhang", min.cells = 0, min.features = 0)
pbmc_zhengzhang@meta.data=pbmc@meta.data[USED_CELL,]
pbmc_zhengzhang <- NormalizeData(pbmc_zhengzhang, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc_zhengzhang)
pbmc_zhengzhang <- ScaleData(pbmc_zhengzhang, features = all.genes)
VariableFeatures(pbmc_zhengzhang)=VariableFeatures(pbmc)
pbmc_zhengzhang <- RunPCA(pbmc_zhengzhang, npcs=10,features = VariableFeatures(object = pbmc_zhengzhang))
pbmc_zhengzhang <- RunUMAP(pbmc_zhengzhang, dims = 1:10)
pbmc_zhengzhang@reductions$umap@cell.embeddings=pbmc@reductions$umap@cell.embeddings[USED_CELL,]
saveRDS(pbmc_zhengzhang,'pbmc_zhengzhang.RDS')




