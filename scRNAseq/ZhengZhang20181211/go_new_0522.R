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




#library(sva)
#library(pamr)
#library(limma)
#pca=pbmc@reductions$pca@cell.embeddings
#batch=as.character(pbmc@meta.data$orig.ident)
#pheno = data.frame(batch=as.matrix(batch))
#edata = t(pca)
#batch = pheno$batch
#modcombat = model.matrix(~1, data=pheno)
#combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
#ttt=t(combat_edata)
#colnames(ttt)=colnames(pbmc@reductions$pca@cell.embeddings)
#rownames(ttt)=rownames(pbmc@reductions$pca@cell.embeddings)
#pbmc@reductions$pca@cell.embeddings=ttt

source('BEER_Seurat3.R')

EXP=as.matrix(pbmc@assays$RNA@counts)
D1=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='CDC42HET'),]) ) ]
D2=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='SmallIntestine'),]) ) ]

one1=.data2one(D1, VariableFeatures(object = pbmc), CPU=4, PCNUM=150, SEED=123,  PP=30)
one2=.data2one(D2, VariableFeatures(object = pbmc), CPU=4, PCNUM=150, SEED=123,  PP=30)
D1X=one1
D2X=one2
CNUM=15
G1=.getGroup(D1X,'D1',CNUM)
G2=.getGroup(D2X,'D2',CNUM)
GROUP=c(G1,G2)
pbmc@meta.data$group=GROUP
VP_OUT=.getValidpair(D1, G1, D2, G2, 4, method='kendall', 10)
VP=VP_OUT$vp
NROW_VP=nrow(VP)
print('n(Validpair):')
print(NROW_VP)
    #if(NROW_VP<=1 | is.null(NROW_VP) ){print('Please set a smaller CNUM !!!')}
    #if(NROW_VP<=1 | is.null(NROW_VP) ){return(message("Please set a smaller CNUM !!!"))}
    ##########################
VPCOR=0
VP=VP[which(VP_OUT$cor>=VPCOR),]
MAP=rep('NA',length(GROUP))
MAP[which(GROUP %in% VP[,1])]='D1'
MAP[which(GROUP %in% VP[,2])]='D2'
pbmc@meta.data$map=MAP

######

BATCH=c(rep('D1',ncol(D1)),rep('D2',ncol(D2)))
DR=pbmc@reductions$pca@cell.embeddings 
B1index=which(BATCH=='D1')
B2index=which(BATCH=='D2')
OUT=.evaluateBatcheffect(DR, B1index, B2index, GROUP, VP)

#############

mybeer=OUT
plot(mybeer$cor, xlab='PCs', ylab="COR", pch=16)
PCUSE <- which(mybeer$cor> 0.65 & mybeer$fdr<0.05 )
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction = "umap")






mybeer=BEER(D1, D2, CNUM=100, PCNUM=300, CPU=4, MTTAG="^mt-", REGBATCH=TRUE)

#par(mfrow=c(1,2))
plot(mybeer$cor, xlab='PCs', ylab="COR", pch=16)
#plot(-log(max(mybeer$fdr,0.00000000001),10), xlab='PCs', ylab='-log10(FDR)', pch=16)


#pbmc <- mybeer$seurat
PCUSE <- which(mybeer$cor> 0.6 & mybeer$fdr<0.05 )
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",group.by='map')







PCUSE=1:150
pbmc <- RunUMAP(pbmc, dims = PCUSE)




DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc,file='pbmc.RDS')






