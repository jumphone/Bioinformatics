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


saveRDS(pbmc,file='pbmc.RDS')

pdf('TYPE.pdf',width=10,height=7)
DimPlot(pbmc, reduction = "umap",label=T)
DimPlot(pbmc, reduction = "umap",group.by='batch',label=T)

dev.off()

pbmc@meta.data$type=paste0(pbmc@meta.data$batch,'_',as.character(pbmc@active.ident))
FeaturePlot(pbmc, features = c("CD3G",'CD3D','CD3E','CD4','FOXP3'))

#show_gene=c('CD3g', 'CD4','FOXP3','CDC42','WAS','GATA3',"CA1","IL4","IL17")

#VlnPlot(pbmc, features = show_gene,slot='data',group.by='type',idents=c('T_Tcell','N_Tcell'))



used=which(pbmc@meta.data$type %in% c('T_Tcell','N_Tcell'))
EXP=as.matrix(pbmc@assays$RNA@data[,used])
colnames(EXP)=paste0(BATCH[used],'_',colnames(EXP))



CD3G=which(rownames(EXP) %in% c('CD3G'))
CD3D=which(rownames(EXP) %in% c('CD3D'))
CD3E=which(rownames(EXP) %in% c('CD3E'))


CD4=which(rownames(EXP) %in% c('CD4'))
FOXP3=which(rownames(EXP) %in% c('FOXP3'))

TMP=EXP[cbind(CD3G,CD3D,CD3E,CD4,FOXP3),]
TMP[which(TMP>0)]=1
#which()

#PAT=EXP[which(rownames(EXP) %in% show_gene),]
heatmap(TMP,scale='none',margins=c(15,10),Colv=F)


table(pbmc@meta.data$type)
#T_Tcell=34
#N_Tcell=11
MIN=apply(TMP[c(4,5),],2,min)

POS=which(MIN>0)
NEG=which(MIN==0)
EXP=cbind(EXP[,NEG],EXP[,POS])
colnames(EXP)[ncol(EXP)]=names(POS)


#EXP[,which(colnames(EXP)==POS)]
PT=rep('NEG',ncol(EXP))
PT[which(colnames(EXP)==names(POS))]='POS'
PT=t(as.matrix(PT))


OUT=cbind(rownames(EXP),rep('NO',nrow(EXP)),EXP)
colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
write.table(OUT,'EXP.txt',sep='\t',quote=F,row.names=F,col.names=T)
write.table(PT,'PT.cls',sep=' ',quote=F,row.names=F,col.names=F )


