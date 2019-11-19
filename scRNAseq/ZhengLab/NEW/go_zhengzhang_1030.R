
setwd('F:/Zhenglab/NewZhengZhang')
library(Seurat)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#REF=readRDS('REF.RDS')
#LABEL=readRDS('LABEL.RDS')



CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO <- Read10X(data.dir = "./CDC42_KO")
CDC42Rescue <- Read10X(data.dir = "./CDC42Rescue")
YapHet <- Read10X(data.dir = "./YapHet")

BATCH=c( rep('CDC42HET',ncol(CDC42HET )), rep('CDC42KO',ncol(CDC42KO )),
        rep('CDC42Rescue',ncol(CDC42Rescue )), rep('YapHet',ncol(YapHet )))


D1=.simple_combine(CDC42HET, CDC42KO)$combine
D2=.simple_combine(CDC42Rescue, YapHet)$combine

DATA=.simple_combine(D1, D2)$combine

############
#QC
#pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
#Idents(pbmc)=BATCH
#pbmc@meta.data$batch=BATCH
#pbmc@meta.data$tag=TAG
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)
#BATCH=pbmc@meta.data$batch
#TAG=pbmc@meta.data$tag
#DATA=as.matrix(pbmc@assays$RNA@counts[,which(colnames(pbmc@assays$RNA@counts) %in% colnames(pbmc@assays$RNA@data))])

#############################

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE )

# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

pbmc <- mybeer$seurat
PCUSE=mybeer$select   
#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)


saveRDS(mybeer,file='mybeer.RDS')
saveRDS(pbmc,file='pbmc.RDS')

##############
DimPlot(pbmc, reduction = "umap", split.by = "batch",ncol=2)

########################




VEC=pbmc@reductions$umap@cell.embeddings

# Here, we use K-means to do the clustering
N=200
set.seed(123)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)+ NoLegend()


pbmc@meta.data$celltype=rep('Enterocyte',ncol(pbmc))
######################

FeaturePlot(pbmc, ncol=3, features=c('Lyz1','Defa17','Ang4','Defa22','Defa24'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('16','48'))]='Paneth.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

##############################

FeaturePlot(pbmc, ncol=3, features=c('Lgr5','Ascl2','Slc12a2','Axin2','Olfm4','Gkn3'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('142','56','116','62'))]='Stem.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

table(pbmc@meta.data$celltype, pbmc@meta.data$batch)

TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)


##############################

FeaturePlot(pbmc, ncol=3, features=c('Mki67','Cdk4','Mcm5','Mcm6','Pcna'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(179,19,27,132,19,20,68,182,24,39,37,51,60,
                                                         67,146,149,
                                                         191,23,148,
                                                         174,104,151,192,178,
                                                         65,172))]='TA.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)

##############################

FeaturePlot(pbmc, ncol=3, features=c('Muc2','Clca3','Tff3','Agr2'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(32,150,126,73,81))]='Goblet.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()




##############################

FeaturePlot(pbmc, ncol=3, features=c('Chga','Chgb','Tac1','Tph1','Neurog3'))


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('183'))]='Endocrine.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()



##############################
FeaturePlot(pbmc, ncol=3, features=c('Dclk1','Trpm5','Gfi1b','Il25'))



pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('123'))]='Tuft.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()





##############################
FeaturePlot(pbmc, ncol=3, features=c('Ptprc','Cd3g'))
pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(49,197,31,167,110,180,34,166,157,
                                                          92,58,36,139,196,141,129,5,103,72,186,145,165,
                                                          15,102,26,80,136,170,107,120,66,94,171,64,91,114,143,22
                                                         ))]='Immune.Cell'


DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()





##############################
FeaturePlot(pbmc, ncol=3, features=c('Alpi','Apoa1','Apoa4','Fabp1'))
pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(166,157,31,110
                                                         ))]='Enterocyte'




DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()


TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)



FeaturePlot(pbmc, ncol=3, features=c('Slc16a1'))






saveRDS(mybeer,file='mybeer.RDS')
saveRDS(pbmc,file='pbmc.RDS')
saveRDS(pbmc,file='pbmc.final.RDS')

############################################################


setwd('F:/Zhenglab/NewZhengZhang')

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS(file='pbmc.final.RDS')


#######################################
ALL.DATA=as.matrix(pbmc@assays$RNA@data)
ALL.TAG=pbmc@meta.data$batch


#######################################
CT='TA.Cell'
BT=c('CDC42HET','CDC42KO')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################

#######################################
CT='Stem.Cell'
BT=c('CDC42HET','CDC42KO')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################




#######################################
CT='TA.Cell'
BT=c('CDC42Rescue','CDC42KO')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################

#######################################
CT='Stem.Cell'
BT=c('CDC42Rescue','CDC42KO')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################


