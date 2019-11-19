
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























########################################
setwd('F:/Zhenglab/NewZhengZhang')
mybeer=readRDS(file='mybeer.RDS')


PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


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

N=20
this_step=1/N

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
            
            if(length(this_in_index)>0){
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

TIME_CUT=2

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


visa.plot3d(VEC.E,COL)

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

