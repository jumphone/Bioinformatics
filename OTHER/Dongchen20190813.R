setwd('./Desktop/CCHMC_Project/Dongchen/DATA')

source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')



DATA <- .readTable(PATH='data.txt', SEP='\t')

.norm_sum=function(x){

    y=x/sum(x)
    return(y)
    }
rownames(DATA)
NDATA=DATA

NDATA[c(8:12),]=apply(DATA[c(8:12),],2,.norm_sum)
NDATA[c(13:17),]=apply(DATA[c(13:17),],2,.norm_sum)
NDATA[c(18:22),]=apply(DATA[c(18:22),],2,.norm_sum)
NDATA[c(24:28),]=apply(DATA[c(24:28),],2,.norm_sum)
NDATA[c(29:33),]=apply(DATA[c(29:33),],2,.norm_sum)




NDATA=t(apply(NDATA ,1,.norm_sum))


library(dplyr)
library(Seurat)
pbmc <- CreateSeuratObject(counts = NDATA, project = "DC", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)




all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


Idents(pbmc)=c(rep('WT',3),rep('HET',5),rep('CKO',9))

VariableFeatures(object = pbmc)=all.genes

#VAR=apply(pbmc@assays$RNA@data,1,var)

#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 8)
VariableFeatures(pbmc)=all.genes
pbmc <- RunPCA(pbmc, npcs =16 ,features =  VariableFeatures(pbmc))


DimPlot(pbmc, reduction = "pca")


ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:10)



DimPlot(pbmc, reduction = "umap")

UMAP1=pbmc@reductions$umap@cell.embeddings[,1]
boxplot(UMAP1~Idents(pbmc))
t.test(UMAP1[which(Idents(pbmc)=='CKO')], UMAP1[which(Idents(pbmc) %in% c('WT','HET'))])



PC1=pbmc@reductions$pca@cell.embeddings[,1]
boxplot(PC1~Idents(pbmc))

t.test(PC1[which(Idents(pbmc)=='CKO')], PC1[which(Idents(pbmc) %in% c('WT','HET'))])
wilcox.test(PC1[which(Idents(pbmc)=='CKO')], PC1[which(Idents(pbmc) %in% c('WT','HET'))])






#DimPlot(pbmc, reduction = "pca")



PC1=pbmc@reductions$pca@cell.embeddings[,1]
boxplot(PC1~Idents(pbmc))

t.test(PC1[which(Idents(pbmc)=='CKO')], PC1[which(Idents(pbmc) %in% c('WT','HET'))])
wilcox.test(PC1[which(Idents(pbmc)=='CKO')], PC1[which(Idents(pbmc) %in% c('WT','HET'))])




pbmc@reductions$pca@feature.loadings[,1]
