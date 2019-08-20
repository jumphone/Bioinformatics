


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


DimPlot(pbmc_zhengzhang, reduction.use='umap', group.by='batch', pt.size=0.1,label=T)




pbmc_zhengzhang=readRDS('pbmc_zhengzhang.RDS')



#pbmc@meta.data$celltype=rep(NA,length(pbmc@meta.data$batch))
#pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='RNA')]=pbmc.rna@meta.data$celltype
#DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)

#######
VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)
N=150
K=kmeans(VEC,centers=N)
pbmc@meta.data$kclust=K$cluster   

pdf('ZZ_CLUST.pdf',width=10,height=10)
DimPlot(pbmc, reduction.use='umap', group.by='kclust', pt.size=0.1,label=T)+NoLegend()
dev.off()
pbmc@meta.data$transfer=rep(NA, length(pbmc@meta.data$celltype))
TMP=cbind(pbmc@meta.data$celltype, pbmc@meta.data$kclust)

KC=unique(pbmc@meta.data$kclust)
i=1
while(i<=length(KC)){
    this_kc=KC[i]
    this_index=which(pbmc@meta.data$kclust==this_kc)
    this_tb=table(pbmc@meta.data$celltype[this_index])
    if(length(this_tb)!=0){
        this_ct=names(this_tb)[which(this_tb==max(this_tb))[1]]
        pbmc@meta.data$transfer[this_index]=this_ct}
    i=i+1}
    
#pbmc@meta.data$tf.ct=pbmc@meta.data$celltype
#NA.index=which(is.na(pbmc@meta.data$celltype))
#pbmc@meta.data$tf.ct[NA.index]=pbmc@meta.data$transfer[NA.index]


DimPlot(pbmc, reduction.use='umap', group.by='transfer', pt.size=0.1,label=T)

################################################################################




pbmc@meta.data$level1=rep('Enterocyte',ncol(pbmc))

pbmc@meta.data$level1[which(pbmc@meta.data$transfer =='Stem')]='Stem'
pbmc@meta.data$level1[which(pbmc@meta.data$transfer %in% c('TA.Early','TA.G1','TA.G2'))]='TA'

pbmc@meta.data$level1[which(pbmc@meta.data$kclust %in% c(145))]='Stem'
pbmc@meta.data$level1[which(pbmc@meta.data$kclust %in% c(36,75,119))]='Endocrine'
pbmc@meta.data$level1[which(pbmc@meta.data$kclust %in% c(141,123,47,138,116))]='Goblet'
pbmc@meta.data$level1[which(pbmc@meta.data$kclust %in% c(105,25))]='Tuft'
pbmc@meta.data$level1[which(pbmc@meta.data$kclust %in% c(103,32,54,128,66,4,99,45,117,115,150,120,18,76,79,191,112,89,121,101))]='Immune'
pbmc@meta.data$level1[which(pbmc@meta.data$kclust %in% c(22,64))]='Paneth'

saveRDS(pbmc@meta.data,'pbmc1_meta.RDS')


#pbmc@meta.data$level1[which(pbmc@meta.data$transfer %in% c('TA.Early','TA.G1','TA.G2'))]='TA'


DimPlot(pbmc, reduction.use='umap', group.by='level1', pt.size=0.1,label=T)



USED_CELL=which(pbmc@meta.data$batch %in% c('CDC42KO','CDC42HET'))
pbmc_zhengzhang@meta.data=pbmc@meta.data[USED_CELL,]
saveRDS(pbmc_zhengzhang@meta.data,'pbmc_zhengzhang_meta.RDS')










