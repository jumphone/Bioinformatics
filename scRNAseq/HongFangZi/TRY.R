
setwd('F:/Hongfangzi')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

DATA=readRDS( file='./DATA.RDS')
BATCH=readRDS( file='./BATCH.RDS')

USED=which(BATCH %in% c('decidua0117','decidua0417.2',
                       'decidua2018c','decidua20190215',
                       'decidua20190420', 'decidua2019c',
                       'decidua508','decidua510'
                       ))

DATA=DATA[,USED]
BATCH=BATCH[USED]

PosN=apply(DATA,2,.getPos)
USED=which(PosN>500 & PosN<2000)

DATA=DATA[,USED]
BATCH=BATCH[USED]


##########################
mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   


# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


saveRDS(mybeer,'mybeer_decidua.RDS')



pbmc=mybeer$seurat
PCUSE <- mybeer$select
#pbmc=BEER.combat(pbmc) 
#umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
#pbmc@reductions$umap@cell.embeddings=umap

pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE,
                n.neighbors = 20, min.dist=0.1)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)

FeaturePlot(pbmc,features='CD34',order=TRUE)
FeaturePlot(pbmc,features='COL1A1',order=TRUE)
FeaturePlot(pbmc,features='VIM',order=TRUE)


saveRDS(pbmc,'pbmc_decidua.RDS')






















FeaturePlot(pbmc,features='CCL5',order=TRUE)


#pbmc=BEER.combat(pbmc)




pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE,
                n.neighbors = 20, min.dist=0.1)


FeaturePlot(pbmc,features='CCL5',order=TRUE)


FeaturePlot(pbmc,features='COL1A1')

FeaturePlot(pbmc,features='KRT7',order=TRUE)

FeaturePlot(pbmc,features='LTB',order=TRUE)

FeaturePlot(pbmc,features='HLA-DPA1',order=TRUE)


FeaturePlot(pbmc,features='CD34')


DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1) 

TYPE=rep('NM',length(pbmc@meta.data$batch))
TYPE[which(pbmc@meta.data$batch %in% c('decidua2019c',
                                      'decidua20190215',
                                      'decidua510'))]='PE'

pbmc@meta.data$type=TYPE


DimPlot(pbmc, reduction.use='umap', group.by='type', pt.size=0.1) 



USED=vector.SeuratSelect(pbmc)
pbmc=vector.SeuratAddMetaByCell(pbmc, USED)

DimPlot(pbmc, reduction.use='umap', group.by='select', pt.size=0.1) 

Idents(pbmc)=pbmc@meta.data$select
MMM=FindMarkers(pbmc, ident.1 = 'YES', only.pos=TRUE,min.pct = 0.25)

head(MMM,n=20)









#PCA = pbmc@reductions$pca@cell.embeddings
#PCA=vector.rankPCA(PCA)
#pbmc@reductions$pca@cell.embeddings=PCA




#umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
#pbmc@reductions$umap@cell.embeddings=umap
#DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)

#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)


FeaturePlot(pbmc,features='CCL5')
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE[1:50], check_duplicates=FALSE,n.neighbors = 5)


DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1) 



































#######################
setwd('F:/Hongfangzi')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

DATA=readRDS( file='./DATA.RDS')
BATCH=readRDS( file='./BATCH.RDS')

USED=which(BATCH %in% c('placenta0423.2','placenta508.1',
                       'placenta508.2','placenta514.2',
                       'placenta2018c', 'placenta2019c',
                       'placenta20190402'
                       ))




DATA=DATA[,USED]
BATCH=BATCH[USED]

library(Seurat)
# DATA: Expression matrix. Rownames are gene names. Colnames are cell names.
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,'pbmc_placenta.RDS')


VEC = pbmc@reductions$umap@cell.embeddings
rownames(VEC) = colnames(pbmc)
PCA = pbmc@reductions$pca@cell.embeddings

Sys.time()
PCA=vector.rankPCA(PCA)
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
Sys.time()

FeaturePlot(pbmc,features='ACTA2')
################################################################################



setwd('F:/Hongfangzi')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

DATA=readRDS( file='./DATA.RDS')
BATCH=readRDS( file='./BATCH.RDS')

USED=which(BATCH %in% c('decidua0117','decidua0417.2',
                       'decidua2018c','decidua20190215',
                       'decidua20190420', 'decidua2019c',
                       'decidua508','decidua510'
                       ))




DATA=DATA[,USED]
BATCH=BATCH[USED]


mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL)   












library(Seurat)
# DATA: Expression matrix. Rownames are gene names. Colnames are cell names.
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc,'pbmc_decidua.RDS')


VEC = pbmc@reductions$umap@cell.embeddings
rownames(VEC) = colnames(pbmc)
PCA = pbmc@reductions$pca@cell.embeddings

Sys.time()
PCA=vector.rankPCA(PCA)
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
Sys.time()

FeaturePlot(pbmc,features='DIO2')




# FB : ACTA2, COL1A1






