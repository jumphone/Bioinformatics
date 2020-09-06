setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

library(Seurat)
pbmc=readRDS('../pbmc.final.RDS')

EXP=pbmc@assays$RNA@counts
BATCH=pbmc@meta.data$batch


wt.pbmc=readRDS('new.wt.pbmc.rds')
ko.pbmc=readRDS('new.ko.pbmc.rds')

BATCH=c(BATCH,rep('NEW.WT',ncol(wt.pbmc)), rep('NEW.KO',ncol(ko.pbmc)))

DATA=.simple_combine(EXP, wt.pbmc@assays$RNA@counts)$combine
DATA=.simple_combine(DATA, ko.pbmc@assays$RNA@counts)$combine


dim(DATA)
length(BATCH)

saveRDS(DATA, 'DATA.RDS')
saveRDS(BATCH, 'BATCH.RDS')

#################################
mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=3000, SEED=1, COMBAT=TRUE )
saveRDS(mybeer,file='mybeer.RDS')



##########################
setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
library(reticulate)
use_python("C:/Users/cchmc/Anaconda3/python.exe",required=T)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

mybeer=readRDS('mybeer.RDS')

# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

#library(reticulate)
#use_python("C:/Users/cchmc/Anaconda3/python.exe",required=T)

pbmc <- mybeer$seurat
PCUSE=mybeer$select   
#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
saveRDS(pbmc,file='pbmc.RDS')
############################################

setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
library(reticulate)
use_python("C:/Users/cchmc/Anaconda3/python.exe",required=T)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

pbmc=readRDS(pbmc,file='pbmc.RDS')

DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
DimPlot(pbmc, reduction = "umap", split.by = "batch",ncol=3)


##########################################
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

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('71','7'))]='Paneth.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

##############################

FeaturePlot(pbmc, ncol=3, features=c('Lgr5','Ascl2','Slc12a2','Axin2','Olfm4','Gkn3'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('182','82','173','39'))]='Stem.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

table(pbmc@meta.data$celltype, pbmc@meta.data$batch)

TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)


##############################


FeaturePlot(pbmc, ncol=3, features=c('Mki67','Cdk4','Mcm5','Mcm6','Pcna'))


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(35,98,96,3,92,198,29,32,179,36,101,104,55,90,59,102,64,21,143,144,
                                                         42,41,168,134))]='TA.Cell'

DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()


TB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)
.norm_sum=function(x){
    return(round(x/sum(x)*100,2))
    }

apply(TB,2,.norm_sum)

##############################



FeaturePlot(pbmc, ncol=3, features=c('Muc2','Clca3','Tff3','Agr2'))

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(27,122,111,188,49))]='Goblet.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()



##############################

FeaturePlot(pbmc, ncol=3, features=c('Chga','Chgb','Tac1','Tph1','Neurog3'))


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(31,57))]='Endocrine.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()



##############################
FeaturePlot(pbmc, ncol=3, features=c('Dclk1','Trpm5','Gfi1b','Il25'))



pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c('43'))]='Tuft.Cell'
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)+ NoLegend()

##############################
FeaturePlot(pbmc, ncol=3, features=c('Ptprc','Cd3g'))
pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% c(72,4,119,25,9,153,54,192,87,132,136,121,94,152,38,110,139,135,169,34,184,105,
                                                          151,177,100,97,159
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

saveRDS(pbmc,file='pbmc.final.RDS')


######################################################

setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS(file='pbmc.final.RDS')


#######################################
ALL.DATA=as.matrix(pbmc@assays$RNA@data)
ALL.TAG=pbmc@meta.data$batch


#######################################
CT='TA.Cell'
BT=c('NEW.WT','NEW.KO')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################

#######################################
CT='Stem.Cell'
BT=c('NEW.WT','NEW.KO')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################





























































#########################################

HQRef=readRDS('HQRef.rds')
ref_tag=cbind(colnames(pbmc),as.character(pbmc@meta.data$celltype))

exp_ref_mat=as.matrix(pbmc@assays$RNA@data)
exp_sc_mat=as.matrix(wt.pbmc@assays$RNA@data)

out=.get_cor(exp_sc_mat, HQRef, method='spearman',CPU=4, print_step=10)
sc_tag=.get_tag_max(out)

ref_vec=pbmc@reductions$umap@cell.embeddings
saveRDS(ref_vec, 'ref_vec.rds')



pbmc=readRDS('../pbmc.final.RDS')
VG=VariableFeatures(pbmc)
saveRDS(VG, file='VG.rds')

exp_ref_mat[which(rownames(exp_ref_mat) %in% VariableFeatures(pbmc)),]

######################################
setwd('F:/Zhenglab/NewZhengZhang/NEW_20200905')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

library(Seurat)
exp_ref_mat=readRDS('exp_ref_mat.rds')
ref_vec=readRDS('ref_vec.rds')
ref_tag=readRDS('ref_tag.rds')
HQRef=readRDS('HQRef.rds')
VG=readRDS('VG.rds')
wt.pbmc=readRDS('new.wt.pbmc.rds')
ko.pbmc=readRDS('new.ko.pbmc.rds')

exp_sc_mat=as.matrix(wt.pbmc@assays$RNA@data)
sc_tag=readRDS('wt.sc_tag.rds')

out =.vec_projection(exp_sc_mat[which(rownames(exp_sc_mat) %in% VG),], sc_tag, exp_ref_mat[which(rownames(exp_ref_mat) %in% VG),], ref_tag, ref_vec, 
        method='spearman', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, CPU=2, print_step=10)







