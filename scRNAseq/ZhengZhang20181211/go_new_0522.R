#Seurat 3.0

library(dplyr)
library(Seurat)



pbmc.data.1 <- Read10X(data.dir = "./CDC42_HET/")
pbmc.data.2 <- Read10X(data.dir = "./Small_Intestine/")

colnames(pbmc.data.1)=paste0('CDC42HET_',colnames(pbmc.data.1))
colnames(pbmc.data.2)=paste0('SmallIntestine_',colnames(pbmc.data.2))

source('scRef.R')
DATA=.simple_combine(pbmc.data.1,pbmc.data.2)$combine

pbmc <- CreateSeuratObject(counts = DATA, project = "Intestine", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

pdf('QC1.pdf',width=15,height=5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 50)
pbmc@meta.data$batch=pbmc@meta.data$orig.ident



pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = c("percent.mt","nCount_RNA","batch"))
PCNUM=310
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=PCNUM)


#################
source('BEER_Seurat3.R')

EXP=as.matrix(pbmc@assays$RNA@counts)
D1=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='CDC42HET'),]) ) ]
D2=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='SmallIntestine'),]) ) ]


mybeer=BEER(D1, D2, CNUM=20, PCNUM=200, GN=5000, CPU=4, MTTAG="^mt-", REGBATCH=FALSE)



plot(mybeer$cor, xlab='PCs', ylab="COR", pch=16)


npbmc <- mybeer$seurat
PCUSE <- which(mybeer$cor> quantile(mybeer$cor,0.2) & mybeer$fdr<0.05 )
npbmc <- RunUMAP(object = npbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
pdf('BEER.pdf',width=10,heigth=8)
DimPlot(npbmc, reduction = "umap")
DimPlot(npbmc, reduction = "umap",group.by='map')
dev.off()

saveRDS(npbmc,file='pbmc.RDS')











