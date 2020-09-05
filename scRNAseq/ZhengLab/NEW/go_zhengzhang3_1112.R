

setwd('F:/Zhenglab/NewZhengZhang3')
library(Seurat)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#REF=readRDS('REF.RDS')
#LABEL=readRDS('LABEL.RDS')



CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO <- Read10X(data.dir = "./CDC42_KO")
CDC42Rescue <- Read10X(data.dir = "./CDC42Rescue")


BATCH=c( rep('CDC42HET',ncol(CDC42HET )), rep('CDC42KO',ncol(CDC42KO )),
        rep('CDC42Rescue',ncol(CDC42Rescue )))


D1=.simple_combine(CDC42HET, CDC42KO)$combine

DATA=.simple_combine(D1, CDC42Rescue)$combine


#############################

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE )

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

##################################################


######################

FeaturePlot(pbmc, ncol=3, features=c('Lyz1','Defa17','Ang4','Defa22','Defa24'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('68'))]='Paneth.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

##############################


FeaturePlot(pbmc, ncol=3, features=c('Lgr5','Ascl2','Slc12a2','Axin2','Olfm4','Gkn3'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(122,77,129,43))]='Stem.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

table(pbmc@meta.data$celltype, pbmc@meta.data$batch)

TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)


##############################


pbmc@meta.data$celltype=rep('Enterocyte',ncol(pbmc))
pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('68'))]='Paneth.Cell'
pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(122,77,129,43))]='Stem.Cell'


##############################

FeaturePlot(pbmc, ncol=3, features=c('Mki67','Cdk4','Mcm5','Mcm6','Pcna'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(91, 37,87,39,78,
                                                          28,180,31,164,167,99,35,70,163,187,56,71,18,52,
                                                          85,97,195,29,45,177,178,199,112,193,11,34,25,55,
                                                          169,140,137
                                                         
                                                         ))]='TA.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)

##############################




FeaturePlot(pbmc, ncol=3, features=c('Muc2','Clca3','Tff3','Agr2'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(3,105,142,171,111,46,155,64))]='Goblet.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()





##############################

FeaturePlot(pbmc, ncol=3, features=c('Chga','Chgb','Tac1','Tph1','Neurog3'))


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('100'))]='Endocrine.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()


TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)


##############################
FeaturePlot(pbmc, ncol=3, features=c('Dclk1','Trpm5','Gfi1b','Il25'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('179'))]='Tuft.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

##############################
FeaturePlot(pbmc, ncol=3, features=c('Ptprc','Cd3g'))


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(7,115,72,76,98,19,49,9,197))]='Immune.Cell'


DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

#############
ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
used.cells <- CellSelector(plot = ppp)
#########################

pbmc@meta.data$celltype[which(colnames(pbmc) %in% used.cells)]='Immune.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

########################

FeaturePlot(pbmc, ncol=3, features=c('Alpi','Apoa1','Apoa4','Fabp1'))



DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()


TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)

################


saveRDS(pbmc,file='pbmc.final.RDS')


############




DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=F)



