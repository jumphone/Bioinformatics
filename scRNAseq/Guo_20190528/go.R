library(dplyr)
library(Seurat)


DATA=read.table('GSE81861_CRC_tumor_all_cells_FPKM.csv',sep=',',header=T,row.names=1)
dim(DATA)

gettype <- function(x){
    y=unlist(strsplit(x, "_"))[3]
    return(y)
}

CN=colnames(DATA)
TYPE=apply(matrix(CN,ncol=1),1,gettype)
table(TYPE)
colnames(DATA)=paste0(TYPE,'_',CN)

getgene <- function(x){
    y=unlist(strsplit(x, "_"))[2]
    return(y)
}

RN=rownames(DATA)
GENE=apply(matrix(RN,ncol=1),1,getgene)
#table(TYPE)
UG=names(which(table(GENE)==1))

rownames(DATA)[which(GENE %in% UG)]=GENE[which(GENE %in% UG)]
DATA=DATA[which(GENE %in% UG),]


TD=DATA




DATA=read.table('GSE81861_CRC_NM_all_cells_FPKM.csv',sep=',',header=T,row.names=1)
dim(DATA)

gettype <- function(x){
    y=unlist(strsplit(x, "_"))[3]
    return(y)
}

CN=colnames(DATA)
TYPE=apply(matrix(CN,ncol=1),1,gettype)
table(TYPE)
colnames(DATA)=paste0(TYPE,'_',CN)

getgene <- function(x){
    y=unlist(strsplit(x, "_"))[2]
    return(y)
}

RN=rownames(DATA)
GENE=apply(matrix(RN,ncol=1),1,getgene)
#table(TYPE)
UG=names(which(table(GENE)==1))

rownames(DATA)[which(GENE %in% UG)]=GENE[which(GENE %in% UG)]
DATA=DATA[which(GENE %in% UG),]


ND=DATA

source('BEER.R')

DATA=.simple_combine(ND,TD)$combine
BATCH=c(rep('N',ncol(ND)),rep('T',ncol(TD)))





pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$batch=BATCH

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress = c("batch",'nCount_RNA'))


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")







