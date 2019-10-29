
#RNA: /home/zhangfeng/disk/project/MATT

setwd('C:/Users/cchmc/Desktop/MATT')

DATA=readRDS('combined.DATA.RDS')

library(Seurat)
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
dim(pbmc)
VAR=apply(pbmc@assays$RNA@data,1,var)

#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 10000)
plot(sort(VAR))

used_feature=names(which(VAR>0.01))
length(used_feature)

VariableFeatures(pbmc)=used_feature
#all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=10)


Idents(pbmc)=colnames(pbmc)

DimPlot(pbmc, dims = c(1, 2), reduction = "pca",pt.size=5, label=T)
DimPlot(pbmc, dims = c(2, 3), reduction = "pca",pt.size=5, label=T)
DimPlot(pbmc, dims = c(3, 4), reduction = "pca",pt.size=5, label=T)

write.table(pbmc@assays$RNA@data,file='NormalizedMatrixATAC.txt',sep='\t',row.names=T,quote=F,col.names=T)

write.table(pbmc@assays$RNA@data[which(rownames(pbmc) %in% used_feature),],
            quote=F,file='NormalizedMatrixATAC_bigVar.txt',sep='\t',row.names=T,col.names=T)

