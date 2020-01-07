
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


#####################################################
setwd('F:/Zhenglab/NewZhengZhang')

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS(file='pbmc.final.RDS')


FeaturePlot(pbmc,features='Cdc42')

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




#######################################
CT='TA.Cell'
BT=c('CDC42Rescue','CDC42HET')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################

#######################################
CT='Stem.Cell'
BT=c('CDC42Rescue','CDC42HET')
USED.CELL=which(pbmc@meta.data$celltype ==CT & pbmc@meta.data$batch %in% BT)
DATA=ALL.DATA[,USED.CELL]
TAG=ALL.TAG[USED.CELL]
PATH=paste0('GSEA/',CT,'.',paste0(BT,collapse  ='.'))
.getGSEAinput(DATA,TAG,PATH)
#######################################


# 2020.01.07



setwd('F:/Zhenglab/NewZhengZhang')

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS(file='pbmc.final.RDS')

#CDC42HET     CDC42KO CDC42Rescue      YapHet 
#       6152        4786        3449        7391 


tiff(paste0("IMG/F1.CDC42HET.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42HET'  )  ]
DimPlot(pbmc,group.by='celltype',cells=USED.CELL)
dev.off()

tiff(paste0("IMG/F1.CDC42KO.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42KO'  )  ]
DimPlot(pbmc,group.by='celltype',cells=USED.CELL)
dev.off()

tiff(paste0("IMG/F1.CDC42Rescue.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42Rescue'  )  ]
DimPlot(pbmc,group.by='celltype',cells=USED.CELL)
dev.off()


###############
library(ggplot2)
tiff(paste0("IMG/F2.Ereg.CDC42HET.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42HET'  )  ]
pbmc@meta.data$Ereg.EXP=pbmc@assays$RNA@data[which(rownames(pbmc)=='Ereg'),]
pbmc@meta.data$Ereg.EXP[which(!pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell'))]=0

p1=FeaturePlot(pbmc,features=c('Ereg.EXP'),cells=USED.CELL,order=TRUE,pt.size=2,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 2))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

dev.off()


tiff(paste0("IMG/F2.Ereg.CDC42KO.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42KO'  )  ]
pbmc@meta.data$Ereg.EXP=pbmc@assays$RNA@data[which(rownames(pbmc)=='Ereg'),]
pbmc@meta.data$Ereg.EXP[which(!pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell'))]=0

p1=FeaturePlot(pbmc,features=c('Ereg.EXP'),cells=USED.CELL,order=TRUE,pt.size=2,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 2))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

dev.off()



tiff(paste0("IMG/F2.Ereg.CDC42Rescue.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42Rescue'  )  ]
pbmc@meta.data$Ereg.EXP=pbmc@assays$RNA@data[which(rownames(pbmc)=='Ereg'),]
pbmc@meta.data$Ereg.EXP[which(!pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell'))]=0

p1=FeaturePlot(pbmc,features=c('Ereg.EXP'),cells=USED.CELL,order=TRUE,pt.size=2,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 2))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
             
dev.off()





a=length(which(pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell') & pbmc@meta.data$batch=='CDC42HET' ))
b=length(which(pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell') & pbmc@meta.data$batch=='CDC42HET' & 
              pbmc@assays$RNA@data[which(rownames(pbmc)=='Ereg'),]>0) )

b/a
846
61
0.07210402


a=length(which(pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell') & pbmc@meta.data$batch=='CDC42KO' ))
b=length(which(pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell') & pbmc@meta.data$batch=='CDC42KO' & 
              pbmc@assays$RNA@data[which(rownames(pbmc)=='Ereg'),]>0) )

b/a
1325
339
0.2558491


a=length(which(pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell') & pbmc@meta.data$batch=='CDC42Rescue' ))
b=length(which(pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell') & pbmc@meta.data$batch=='CDC42Rescue' & 
              pbmc@assays$RNA@data[which(rownames(pbmc)=='Ereg'),]>0) )

b/a
766
85
0.1109661

             
             
             
             
######################################
YAP=as.character(read.table('IMG/YAP.txt',header=F)[,1])
MTOR=as.character(read.table('IMG/MTOR.txt',header=F)[,1])
PRO=as.character(read.table('IMG/PRO.txt',header=F)[,1])
###############          
             
             
             
             
###########             
THIS='Stem.Cell'
THIS_INDEX_HET=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42HET')
THIS_INDEX_KO=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42KO')
THIS_INDEX_RES=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42Rescue' )

MAT=pbmc@assays$RNA@data[which(rownames(pbmc) %in% VariableFeatures(pbmc)),c(THIS_INDEX_HET,THIS_INDEX_KO,THIS_INDEX_RES)]

             
         
MAT.BATCH=c(rep('Het',length(THIS_INDEX_HET)),
            rep('KO',length(THIS_INDEX_KO)),
            rep('Rescue',length(THIS_INDEX_RES)))


             
############################################################             
             
GENE=toupper(rownames(MAT))           

YAP.INDEX=which(GENE %in% YAP)
MTOR.INDEX=which(GENE %in% MTOR)
PRO.INDEX=which(GENE %in% PRO)
GENE.BATCH=   c(rep('YAP',length(YAP.INDEX)),
            rep('MTOR',length(MTOR.INDEX)),
            rep('PRO',length(PRO.INDEX)))     

MAT=MAT[c(YAP.INDEX,MTOR.INDEX,PRO.INDEX),]
YAP.INDEX=which(GENE.BATCH=='YAP') 
MTOR.INDEX=which(GENE.BATCH=='MTOR')
PRO.INDEX=which(GENE.BATCH=='PRO')
             
########################################################            
             
             
library('ComplexHeatmap')
library('circlize')
library('seriation')
set.seed(1)


.normX <- function(x){
    y=(x-min(x))/(max(x)-min(x))
    return(y)
    }


.normR <- function(x){
    y=rank(x)/max(rank(x,ties.method='average'))
    return(y)
    }

MAT=as.matrix(MAT)   
s.mat=t(apply(t(MAT),2,.normR))

SSS=function(x,d){return(smooth.spline(x,df=d)$y)}
s.mat=t(apply(t(s.mat),2,SSS,d=30))             
             
LW1=apply(s.mat[,which(MAT.BATCH %in% c('Het'))] ,1,mean)
LW2=apply(s.mat[,which(MAT.BATCH %in% c('Rescue'))] ,1,mean) 
UP= apply(s.mat[,which(MAT.BATCH %in% c('KO'))] ,1,mean)  
             
###################################################
SCORE=apply(cbind(LW1/UP,LW2/UP),1,max)
###################################################
USED=c(YAP.INDEX[order(SCORE[YAP.INDEX])[1:10]],
       MTOR.INDEX[order(SCORE[MTOR.INDEX])[1:7]],
       PRO.INDEX[order(SCORE[PRO.INDEX])[1:20]])
             
###################################################             
                        
o.mat=s.mat[USED,]
o.mat.GENE.BATCH=GENE.BATCH[USED]
#o.mat=MAT[1:100,]
col_fun =colorRamp2(c(0,0.25,0.5,0.75,1 ), c('blue3','blue1','white','red1','red3'))

ha = HeatmapAnnotation(
	 
    Batch = MAT.BATCH,
    col = list(
	       Batch = c("Het"='gold','KO'='royalblue3','Rescue'='grey70')
	      )
    #gp = gpar(col = "black")    
)
             
ha.row=rowAnnotation(
    Pathway = o.mat.GENE.BATCH,
    col = list(
	       Pathway = c("YAP"='red','MTOR'='blue','PRO'='green')
	      )
)
        
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)

             


tiff(paste0("IMG/F3.Stem.heat.tiff"),width=8,height=7,units='in',res=600)             
     
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)
dev.off()
###########################################################################################################



             
             
             
             
             
             
             
             
             
###########             
THIS='TA.Cell'
THIS_INDEX_HET=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42HET')
THIS_INDEX_KO=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42KO')
THIS_INDEX_RES=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42Rescue' )

MAT=pbmc@assays$RNA@data[which(rownames(pbmc) %in% VariableFeatures(pbmc)),c(THIS_INDEX_HET,THIS_INDEX_KO,THIS_INDEX_RES)]

             
         
MAT.BATCH=c(rep('Het',length(THIS_INDEX_HET)),
            rep('KO',length(THIS_INDEX_KO)),
            rep('Rescue',length(THIS_INDEX_RES)))


             
############################################################             
             
GENE=toupper(rownames(MAT))           

YAP.INDEX=which(GENE %in% YAP)
MTOR.INDEX=which(GENE %in% MTOR)
PRO.INDEX=which(GENE %in% PRO)
GENE.BATCH=   c(rep('YAP',length(YAP.INDEX)),
            rep('MTOR',length(MTOR.INDEX)),
            rep('PRO',length(PRO.INDEX)))     

MAT=MAT[c(YAP.INDEX,MTOR.INDEX,PRO.INDEX),]
YAP.INDEX=which(GENE.BATCH=='YAP') 
MTOR.INDEX=which(GENE.BATCH=='MTOR')
PRO.INDEX=which(GENE.BATCH=='PRO')
             
########################################################            
             
             
library('ComplexHeatmap')
library('circlize')
library('seriation')
set.seed(1)


.normX <- function(x){
    y=(x-min(x))/(max(x)-min(x))
    return(y)
    }


.normR <- function(x){
    y=rank(x)/max(rank(x,ties.method='average'))
    return(y)
    }

MAT=as.matrix(MAT)   
s.mat=t(apply(t(MAT),2,.normR))

SSS=function(x,d){return(smooth.spline(x,df=d)$y)}
s.mat=t(apply(t(s.mat),2,SSS,d=30))             
             
LW1=apply(s.mat[,which(MAT.BATCH %in% c('Het'))] ,1,mean)
LW2=apply(s.mat[,which(MAT.BATCH %in% c('Rescue'))] ,1,mean) 
UP= apply(s.mat[,which(MAT.BATCH %in% c('KO'))] ,1,mean)  
             
###################################################
SCORE=apply(cbind(LW1/UP,LW2/UP),1,max)
###################################################
USED=c(YAP.INDEX[order(SCORE[YAP.INDEX])[1:10]],
       MTOR.INDEX[order(SCORE[MTOR.INDEX])[1:7]],
       PRO.INDEX[order(SCORE[PRO.INDEX])[1:20]])
             
###################################################             
                        
o.mat=s.mat[USED,]
o.mat.GENE.BATCH=GENE.BATCH[USED]
#o.mat=MAT[1:100,]
col_fun =colorRamp2(c(0,0.25,0.5,0.75,1 ), c('blue3','blue1','white','red1','red3'))

ha = HeatmapAnnotation(
	 
    Batch = MAT.BATCH,
    col = list(
	       Batch = c("Het"='gold','KO'='royalblue3','Rescue'='grey70')
	      )
    #gp = gpar(col = "black")    
)
             
ha.row=rowAnnotation(
    Pathway = o.mat.GENE.BATCH,
    col = list(
	       Pathway = c("YAP"='red','MTOR'='blue','PRO'='green')
	      )
)
        
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)

             


tiff(paste0("IMG/F3.TA.heat.tiff"),width=8,height=7,units='in',res=600)             
     
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)
dev.off()
###########################################################################################################


        
             
    
            
             
             
             
             
             
#AD. AP             
AD=as.character(read.table('IMG/AD.txt',header=F)[,1])
AP=as.character(read.table('IMG/AP.txt',header=F)[,1])            
###########             
THIS='Stem.Cell'
THIS_INDEX_HET=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42HET')
THIS_INDEX_KO=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42KO')
THIS_INDEX_RES=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42Rescue' )

MAT=pbmc@assays$RNA@data[which(rownames(pbmc) %in% VariableFeatures(pbmc)),c(THIS_INDEX_HET,THIS_INDEX_KO,THIS_INDEX_RES)]

             
         
MAT.BATCH=c(rep('Het',length(THIS_INDEX_HET)),
            rep('KO',length(THIS_INDEX_KO)),
            rep('Rescue',length(THIS_INDEX_RES)))


             
############################################################             
             
GENE=toupper(rownames(MAT))           

AD.INDEX=which(GENE %in% AD)
AP.INDEX=which(GENE %in% AP)

GENE.BATCH=   c(rep('AD',length(AD.INDEX)),
            rep('AP',length(AP.INDEX)))     

MAT=MAT[c(AD.INDEX,AP.INDEX),]
AD.INDEX=which(GENE.BATCH=='AD') 
AP.INDEX=which(GENE.BATCH=='AP')
          
########################################################            
             
             
library('ComplexHeatmap')
library('circlize')
library('seriation')
set.seed(1)


.normX <- function(x){
    y=(x-min(x))/(max(x)-min(x))
    return(y)
    }


.normR <- function(x){
    y=rank(x)/max(rank(x,ties.method='average'))
    return(y)
    }

MAT=as.matrix(MAT)   
s.mat=t(apply(t(MAT),2,.normR))

SSS=function(x,d){return(smooth.spline(x,df=d)$y)}
s.mat=t(apply(t(s.mat),2,SSS,d=30))             
             
LW1=apply(s.mat[,which(MAT.BATCH %in% c('Het'))] ,1,mean)
LW2=apply(s.mat[,which(MAT.BATCH %in% c('Rescue'))] ,1,mean) 
UP= apply(s.mat[,which(MAT.BATCH %in% c('KO'))] ,1,mean)  
             
###################################################
#SCORE=apply(cbind(LW1/UP,LW2/UP),1,max)
SCORE=LW1/UP
###################################################
USED=c(AD.INDEX[order(SCORE[AD.INDEX])[1:9]],
       AP.INDEX[order(SCORE[AP.INDEX])[1:8]])
             
###################################################             
                        
o.mat=s.mat[USED,]
o.mat.GENE.BATCH=GENE.BATCH[USED]
#o.mat=MAT[1:100,]
col_fun =colorRamp2(c(0,0.25,0.5,0.75,1 ), c('blue3','blue1','white','red1','red3'))

ha = HeatmapAnnotation(
	 
    Batch = MAT.BATCH,
    col = list(
	       Batch = c("Het"='gold','KO'='royalblue3','Rescue'='grey70')
	      )
    #gp = gpar(col = "black")    
)
             
ha.row=rowAnnotation(
    Pathway = o.mat.GENE.BATCH,
    col = list(
	       Pathway = c("AD"='red','AP'='blue')
	      )
)
        
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)

             


tiff(paste0("IMG/F4.Stem.heat.ADAP.tiff"),width=8,height=7,units='in',res=600)             
     
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)
dev.off()
#####################             
             
             
             
     
             
             
#AD. AP             
AD=as.character(read.table('IMG/AD.txt',header=F)[,1])
AP=as.character(read.table('IMG/AP.txt',header=F)[,1])            
###########             
THIS='TA.Cell'
THIS_INDEX_HET=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42HET')
THIS_INDEX_KO=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42KO')
THIS_INDEX_RES=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42Rescue' )

MAT=pbmc@assays$RNA@data[which(rownames(pbmc) %in% VariableFeatures(pbmc)),c(THIS_INDEX_HET,THIS_INDEX_KO,THIS_INDEX_RES)]

             
         
MAT.BATCH=c(rep('Het',length(THIS_INDEX_HET)),
            rep('KO',length(THIS_INDEX_KO)),
            rep('Rescue',length(THIS_INDEX_RES)))


             
############################################################             
             
GENE=toupper(rownames(MAT))           

AD.INDEX=which(GENE %in% AD)
AP.INDEX=which(GENE %in% AP)

GENE.BATCH=   c(rep('AD',length(AD.INDEX)),
            rep('AP',length(AP.INDEX)))     

MAT=MAT[c(AD.INDEX,AP.INDEX),]
AD.INDEX=which(GENE.BATCH=='AD') 
AP.INDEX=which(GENE.BATCH=='AP')
          
########################################################            
             
             
library('ComplexHeatmap')
library('circlize')
library('seriation')
set.seed(1)


.normX <- function(x){
    y=(x-min(x))/(max(x)-min(x))
    return(y)
    }


.normR <- function(x){
    y=rank(x)/max(rank(x,ties.method='average'))
    return(y)
    }

MAT=as.matrix(MAT)   
s.mat=t(apply(t(MAT),2,.normR))

SSS=function(x,d){return(smooth.spline(x,df=d)$y)}
s.mat=t(apply(t(s.mat),2,SSS,d=30))             
             
LW1=apply(s.mat[,which(MAT.BATCH %in% c('Het'))] ,1,mean)
LW2=apply(s.mat[,which(MAT.BATCH %in% c('Rescue'))] ,1,mean) 
UP= apply(s.mat[,which(MAT.BATCH %in% c('KO'))] ,1,mean)  
             
###################################################
#SCORE=apply(cbind(LW1/UP,LW2/UP),1,max)
SCORE=LW1/UP
###################################################
USED=c(AD.INDEX[order(SCORE[AD.INDEX])[1:14]],
       AP.INDEX[order(SCORE[AP.INDEX])[1:8]])
             
###################################################             
                        
o.mat=s.mat[USED,]
o.mat.GENE.BATCH=GENE.BATCH[USED]
#o.mat=MAT[1:100,]
col_fun =colorRamp2(c(0,0.25,0.5,0.75,1 ), c('blue3','blue1','white','red1','red3'))

ha = HeatmapAnnotation(
	 
    Batch = MAT.BATCH,
    col = list(
	       Batch = c("Het"='gold','KO'='royalblue3','Rescue'='grey70')
	      )
    #gp = gpar(col = "black")    
)
             
ha.row=rowAnnotation(
    Pathway = o.mat.GENE.BATCH,
    col = list(
	       Pathway = c("AD"='red','AP'='blue')
	      )
)
        
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)

             


tiff(paste0("IMG/F4.TA.heat.ADAP.tiff"),width=8,height=7,units='in',res=600)             
     
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)
dev.off()
#####################             
                         
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
###########             
THIS='TA.Cell'
THIS_INDEX_HET=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42HET')
THIS_INDEX_KO=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42KO')
THIS_INDEX_RES=which(pbmc@meta.data$celltype %in% c(THIS) & pbmc@meta.data$batch=='CDC42Rescue' )

MAT=pbmc@assays$RNA@data[which(rownames(pbmc) %in% VariableFeatures(pbmc)),c(THIS_INDEX_HET,THIS_INDEX_KO,THIS_INDEX_RES)]

             
         
MAT.BATCH=c(rep('Het',length(THIS_INDEX_HET)),
            rep('KO',length(THIS_INDEX_KO)),
            rep('Rescue',length(THIS_INDEX_RES)))

           
GENE=toupper(rownames(MAT))           

YAP.INDEX=which(GENE %in% YAP)
MTOR.INDEX=which(GENE %in% MTOR)
PRO.INDEX=which(GENE %in% PRO)
MAT=MAT[c(YAP.INDEX,MTOR.INDEX,PRO.INDEX),]


.normR <- function(x){
    y=rank(x)/max(rank(x,ties.method='average'))
    return(y)
    }

MAT=as.matrix(MAT)   
s.mat=t(apply(t(MAT),2,.normR))

SSS=function(x,d){return(smooth.spline(x,df=d)$y)}
s.mat=t(apply(t(s.mat),2,SSS,d=30))          
o.mat=s.mat[USED,]         
             
o.mat.GENE.BATCH=GENE.BATCH[USED]
#o.mat=MAT[1:100,]
col_fun =colorRamp2(c(0,0.25,0.5,0.75,1 ), c('blue3','blue1','white','red1','red3'))

ha = HeatmapAnnotation(
	 
    Batch = MAT.BATCH,
    col = list(
	       Batch = c("Het"='gold','KO'='royalblue3','Rescue'='grey70')
	      )
    #gp = gpar(col = "black")    
)
             
ha.row=rowAnnotation(
    Pathway = o.mat.GENE.BATCH,
    col = list(
	       Pathway = c("YAP"='red','MTOR'='blue','PRO'='green')
	      )
)
        
Heatmap(o.mat,row_title='',name="Exp",cluster_rows=FALSE,cluster_columns=FALSE,
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=TRUE,
	col=col_fun, border = TRUE,
	top_annotation = ha,right_annotation=ha.row
	
	)      
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             







########################################
setwd('F:/Zhenglab/NewZhengZhang')
mybeer=readRDS(file='mybeer.RDS')


PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

pbmc <- mybeer$seurat
PCUSE=mybeer$select   
#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10,DM=3)


source('https://raw.githubusercontent.com/jumphone/VISA/master/VISA.R')

COL=visa.col(pbmc@meta.data$batch)

visa.plot3d(umap, COL)




DATA=as.matrix(pbmc@assays$RNA@data)

colnames(DATA)=paste0('C.',1:ncol(DATA))




################################################

VEC=umap
rownames(VEC)=colnames(DATA)

.normOne = function(x){
    x=x
    delta=0.001
    if(var(x)==0){
        y=rep(0,length(x))
    }else{
        y=(x-min(x))/((max(x)-min(x))*(1+delta))
        }
    return(y)
    }

VEC=VEC
VEC.E=apply(VEC,2,.normOne)

visa.plot3d(VEC.E,COL)
####################################################




visa.plot3d(VEC.E,COL)

N=10
this_step=1/N

NUM_CUT=1


SHOW=FALSE

INDEX_LIST=list()
CENTER_LIST=list()

this_x=0
while(this_x<1){
    this_y=0
    while(this_y<1){
        this_z=0
        while(this_z<1){
           
            
            this_in_index = which(VEC.E[,1]>=this_x &  VEC.E[,1] <this_x+this_step &
                                  VEC.E[,2]>=this_y &  VEC.E[,2] <this_y+this_step &
                                  VEC.E[,3]>=this_z &  VEC.E[,3] <this_z+this_step)
            
                
                
            this_center=c(this_x+this_step/2,this_y+this_step/2,this_z+this_step/2)      
            
            if(length(this_in_index)>=NUM_CUT){
                #####################
                INDEX_LIST=c(INDEX_LIST, list(this_in_index))
                CENTER_LIST=c(CENTER_LIST, list(this_center))
                ######################
                if(SHOW==TRUE){
                    points3d(this_center[1],this_center[2],this_center[3],color= 'black')
                    }
                ############
                }
            
            
            this_z=this_z+this_step}            
        this_y=this_y+this_step}
    print(this_x)
    this_x=this_x+this_step}



#############################






DATA=DATA

visa.plot3d(VEC.E,COL)

VEC.E=VEC.E
USED_INDEX=c()
SHOW=TRUE

i=1
while(i<=length(CENTER_LIST)){
    this_center=CENTER_LIST[[i]]
    this_index=INDEX_LIST[[i]]
    ###################
        
    if(length(this_index)>1){    
        this_dist=as.matrix(dist(VEC.E[this_index,]))
        this_sum=apply(this_dist,2,sum)
        ######################
        this_used=names(which(this_sum==min(this_sum))[1])
        this_used_index=which(colnames(DATA)==this_used)
     }else{
        this_used_index=this_index   
        }
    ######################
    if(SHOW==TRUE){
        points3d(VEC.E[this_used_index,1],VEC.E[this_used_index,2],VEC.E[this_used_index,3],
              size=5,color= 'blue')
        }
    ############
    
    USED_INDEX=c(USED_INDEX, this_used_index)
    print(length(USED_INDEX))
    print(i)
    i=i+1}


USED_MAT=DATA[,USED_INDEX]

##################################################



#install.packages('matlib')
library(matlib)

visa.plot3d(VEC.E,COL)



DATA=DATA
CNUM=length(USED_INDEX)
CNAME=colnames(DATA)[USED_INDEX]

TIME_CUT=1.2

SHOW=TRUE

p1=c()
p2=c()
edge_score=c()
EDGE_MAT=c()
L_CUT=5

i=1
while(i<CNUM){
    this_p1_loc=CENTER_LIST[[i]]    
    this_distance_list=c()
    j=i+1
    while(j<=CNUM){
        this_p2_loc=CENTER_LIST[[j]]
        this_distance=sqrt(sum((this_p1_loc-this_p2_loc)^2))
        this_distance_list=c(this_distance_list, this_distance)
        j=j+1}
    
    used_j=i+which(this_distance_list<TIME_CUT*this_step)
    for(j in used_j){
        this_p2_loc=CENTER_LIST[[j]]
        this_p1=CNAME[i]
        this_p2=CNAME[j]
        #################################
        ###################################    
        index_i=INDEX_LIST[[i]]
        index_j=INDEX_LIST[[j]]
        
        if(length(c(index_i,index_j)) >=L_CUT){
      
        direction_ij=this_p1_loc-this_p2_loc

        vec_ij=VEC.E[c(index_i,index_j),]    
        data_ij= DATA[,c(index_i,index_j)]   
         
            
        library(matlib)
        p_list=c()
        iii=1
        while(iii<=nrow(vec_ij)){
            this_p=sqrt(sum(vec_ij[iii,]^2))*cos(angle(as.vector(direction_ij), as.vector( vec_ij[iii,] ),degree=FALSE))
            p_list=c(p_list, this_p)
            iii=iii+1}
        
        this_var=apply(data_ij,1,var)   
        cor_list=rep(0,length(this_var))
        cor_list    
        this_corf=function(x){
            return(cor(x,p_list,method='spearman') )
            }
        cor_list[which(this_var>0)]=apply(data_ij[which(this_var>0),],1,this_corf) 
        names(cor_list)=rownames(data_ij)   
            
            
        ####################
        ########################
        this_score=sqrt(sum((this_p1_loc-this_p2_loc)^2))
        ####################
        if(SHOW==TRUE){
            segments3d(c(this_p1_loc[1],this_p2_loc[1]),
                       c(this_p1_loc[2],this_p2_loc[2]),
                       c(this_p1_loc[3],this_p2_loc[3]),
                       col='black')
            }
        ######################
        p1=c(p1,this_p1)
        p2=c(p2,this_p2) 
        edge_score=c(edge_score, this_score)
        EDGE_MAT=cbind(EDGE_MAT, cor_list)
        ###################
        ################
        }   
        ############
        }
        
        
        
        print(i)
        i=i+1
    }


sort(-apply(abs(EDGE_MAT),1,sum))[1:10]





CUT=0.3
CSIZE=c()
i=1
while(i <= nrow(EDGE_MAT)){



    NET=cbind(p1,p2)
    NET=NET[which(abs(EXP)>=CUT),]
    g <- make_graph(t(NET), directed = FALSE)

    EXP=EDGE_MAT[i,]
    CSIZE=c(CSIZE,max(components(g)$csize))
    if(i %%100==1){print(i)}
    i=i+1
    }


names(CSIZE)=rownames(EDGE_MAT)

    NET=cbind(p1,p2)
    NET=NET[which(abs(EXP)>=CUT),]
    g <- make_graph(t(NET), directed = FALSE)


visa.plot3d(VEC.E,COL)

A=apply(abs(EDGE_MAT),2,sum)

j=1
while(j<=100){

    MAXA=which(A== -sort(-A)[j])

    this_used=which(rownames(VEC.E) %in% c(p1[MAXA],p2[MAXA]))

    points3d(VEC.E[this_used,1],VEC.E[this_used,2],VEC.E[this_used,3],
              size=10,color= 'red')

    j=j+1
}


.writeTable(OUT,PATH='OK.txt')

p1
p2

library(igraph)
NET = cbind(p1,p2) 
g <- make_graph(t(NET), directed = FALSE)
#MST=mst(g, weights = edge_score, algorithm = NULL)























































##################################


visa.plot3d(VEC.E,COL)

VEC.E=VEC.E
USED_INDEX=c()
SHOW=TRUE

i=1
while(i<=length(CENTER_LIST)){
    this_center=CENTER_LIST[[i]]
    this_index=INDEX_LIST[[i]]
    ###################
        
    if(length(this_index)>1){    
        this_dist=as.matrix(dist(VEC.E[this_index,]))
        this_sum=apply(this_dist,2,sum)
        ######################
        this_used=names(which(this_sum==min(this_sum))[1])
        this_used_index=which(colnames(DATA)==this_used)
     }else{
        this_used_index=this_index   
        }
    ######################
    if(SHOW==TRUE){
        points3d(VEC.E[this_used_index,1],VEC.E[this_used_index,2],VEC.E[this_used_index,3],
              size=5,color= 'blue')
        }
    ############
    
    USED_INDEX=c(USED_INDEX, this_used_index)
    print(length(USED_INDEX))
    print(i)
    i=i+1}


USED_MAT=DATA[,USED_INDEX]

##################################################



visa.plot3d(VEC.E,COL)

CNUM=length(USED_INDEX)
CNAME=colnames(DATA)[USED_INDEX]

TIME_CUT=1.2

SHOW=TRUE

p1=c()
p2=c()
edge_score=c()
i=1
while(i<CNUM){
    this_p1_loc=CENTER_LIST[[i]]    
    this_distance_list=c()
    j=i+1
    while(j<=CNUM){
        this_p2_loc=CENTER_LIST[[j]]
        this_distance=sqrt(sum((this_p1_loc-this_p2_loc)^2))
        this_distance_list=c(this_distance_list, this_distance)
        j=j+1}
    
    used_j=i+which(this_distance_list<TIME_CUT*this_step)
    for(j in used_j){
        this_p2_loc=CENTER_LIST[[j]]
        this_p1=CNAME[i]
        this_p2=CNAME[j]
        
        ####################
        #this_cor=cor(USED_MAT[,i],USED_MAT[,j],method='spearman')
        #this_score=1-this_cor 
        ########################
        this_score=sqrt(sum((this_p1_loc-this_p2_loc)^2))
        ####################
        if(SHOW==TRUE){
            segments3d(c(this_p1_loc[1],this_p2_loc[1]),
                       c(this_p1_loc[2],this_p2_loc[2]),
                       c(this_p1_loc[3],this_p2_loc[3]),
                       col='black')
            }
        ######################
        p1=c(p1,this_p1)
        p2=c(p2,this_p2) 
        edge_score=c(edge_score, this_score)
        ###################
        }
        print(i)
        i=i+1
    }

###################################

USED.VEC=VEC.E[USED_INDEX,]

######################

EXP=DATA[which(rownames(DATA)=='Lgr5'),]
EXP=EXP/max(EXP)
COL.E=visa.vcol(EXP, c(0,0.5,1), c('grey90','red1','red3'))
###########
visa.plot3d(VEC.E,COL.E)
#visa.id3d(VEC.E)
######################





library(igraph)
NET = cbind(p1,p2) 
g <- make_graph(t(NET), directed = FALSE)
MST=mst(g, weights = edge_score, algorithm = NULL)

MST.EL=as_edgelist(MST)

SHOW=TRUE

i=1
while(i<=nrow(MST.EL)){
    this_p1=MST.EL[i,1]
    this_p2=MST.EL[i,2]
    this_p1_index=which(CNAME==this_p1)
    this_p2_index=which(CNAME==this_p2)
    #this_p1_loc=CENTER_LIST[[this_p1_index]]
    #this_p2_loc=CENTER_LIST[[this_p2_index]]
    this_p1_loc=USED.VEC[this_p1_index,]
    this_p2_loc=USED.VEC[this_p2_index,]
    ####################
    if(SHOW==TRUE){
        segments3d(c(this_p1_loc[1],this_p2_loc[1]),
                   c(this_p1_loc[2],this_p2_loc[2]),
                   c(this_p1_loc[3],this_p2_loc[3]),
                   col='purple')
        }
    ######################    
    i=i+1}






####################


BTN=betweenness(MST, v = V(MST), directed = TRUE, weights = NULL,
  nobigint = TRUE, normalized = FALSE)
plot(sort(BTN),type='s')
#####################################
BTN.CUT=quantile(BTN,0.5)
#BTN.CUT=1000
USED.NODE=names(which(BTN>=BTN.CUT))
##################


COL.G=rep('grey90',nrow(VEC.E))
visa.plot3d(VEC.E,COL.G)
MST.EL.USED=MST.EL[which(MST.EL[,1] %in% USED.NODE & MST.EL[,2] %in% USED.NODE),]


#############################
SHOW=TRUE

i=1
while(i<=nrow(MST.EL.USED)){
    this_p1=MST.EL.USED[i,1]
    this_p2=MST.EL.USED[i,2]
    this_p1_index=which(CNAME==this_p1)
    this_p2_index=which(CNAME==this_p2)
    #this_p1_loc=CENTER_LIST[[this_p1_index]]
    #this_p2_loc=CENTER_LIST[[this_p2_index]]
    this_p1_loc=USED.VEC[this_p1_index,]
    this_p2_loc=USED.VEC[this_p2_index,]
    ####################
    if(SHOW==TRUE){
        segments3d(c(this_p1_loc[1],this_p2_loc[1]),
                   c(this_p1_loc[2],this_p2_loc[2]),
                   c(this_p1_loc[3],this_p2_loc[3]),
                   col='purple')
        }
    ######################    
    i=i+1}
##################################################

points3d(VEC.E[1,1],VEC.E[1,2],VEC.E[1,3],size=10,color='red')

