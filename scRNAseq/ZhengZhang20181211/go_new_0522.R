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
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress = c("percent.mt","nCount_RNA","batch"))

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))



