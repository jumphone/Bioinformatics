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


#################
source('BEER_Seurat3.R')

EXP=as.matrix(pbmc@assays$RNA@counts)
D1=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='CDC42HET'),]) ) ]
D2=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='SmallIntestine'),]) ) ]


mybeer=BEER(D1, D2, CNUM=20, PCNUM=200, GN=5000, CPU=4, MTTAG="^mt-", REGBATCH=FALSE)



plot(mybeer$cor, xlab='PCs', ylab="COR", pch=16)


npbmc <- mybeer$seurat
PCUSE <- which(mybeer$cor> quantile(mybeer$cor,0.2) & mybeer$fdr<0.05 )
npbmc <- RunUMAP(object = npbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
pdf('BEER.pdf',width=10,height=8)
DimPlot(npbmc, reduction = "umap")
DimPlot(npbmc, reduction = "umap",group.by='map')
dev.off()

saveRDS(npbmc,file='pbmc.RDS')

#############################
VEC=npbmc@reductions$umap@cell.embeddings
META=npbmc@meta.data
saveRDS(VEC, file=paste0('./VEC.RDS'))
saveRDS(META, file=paste0('./META.RDS'))

#########################################
##########################################




##################

VEC=readRDS('VEC.RDS')


NUMBER=150

x_step=(max(VEC[,1])-min(VEC[,1]))/NUMBER
y_step=(max(VEC[,2])-min(VEC[,2]))/NUMBER

X=VEC[,1]
Y=VEC[,2]


############################
tiff('GRID1.tiff',width=9,height=9, res=400,units='in')
plot(VEC,pch=16,cex=0.3,col='grey30')

i=min(X)
while(i<=max(X)){
    abline(v=i,col='grey70')  
    i=i+x_step}

i=min(Y)
while(i<=max(Y)){
    abline(h=i,col='grey70')
    i=i+y_step}
dev.off()
#######################


CELL_X=rep(NA,nrow(VEC))
CELL_Y=rep(NA,nrow(VEC))
CELL_VEC=c()
CELL_X_INDEX=c()
CELL_Y_INDEX=c()

i=0
vec_x=min(X)+x_step*i

while(vec_x<=max(X)){
    j=0
    vec_y=min(Y)+y_step*j

    while(vec_y<=max(Y)){ 
      
        lbx=vec_x
        lby=vec_y
        rtx=vec_x+x_step
        rty=vec_y+y_step
        #rect(lbx,lby,rtx,rty)
      
        CELL_VEC=cbind(CELL_VEC, c((lbx+rtx)/2 , (lby+rty)/2))
        CELL_X_INDEX=c(CELL_X_INDEX,i)
        CELL_Y_INDEX=c(CELL_Y_INDEX,j)
        
        this_cell=which(X>=lbx & X<rtx & Y>=lby & Y<rty)
        if(length(this_cell)>0){        
            CELL_X[this_cell]=i
            CELL_Y[this_cell]=j
        }
             
        j=j+1
        vec_y=min(Y)+y_step*j

        }
    print(i)
    i=i+1
    vec_x=min(X)+x_step*i
}




CELL_VEC=t(CELL_VEC)
rownames(CELL_VEC)=paste0('X_',as.character(CELL_X_INDEX),'_Y_',as.character(CELL_Y_INDEX))
colnames(CELL_VEC)=c('UMAP_1','UMAP_2')


CELL_LABEL=paste0('X_',as.character(CELL_X),'_Y_',as.character(CELL_Y))
V_CELL_VEC=CELL_VEC[which(rownames(CELL_VEC) %in% CELL_LABEL),]

tiff('GRID2.tiff',width=9,height=9, res=400,units='in')
plot(V_CELL_VEC,pch=15,col='grey50',cex=0.4)
dev.off()



NUM=c()
i=1
while(i<=nrow(V_CELL_VEC)){
    #this_vec=V_CELL_VEC[i,]
    this_name=rownames(V_CELL_VEC)[i]
    this_num=length(which(CELL_LABEL ==this_name))
    NUM=c(NUM,this_num)   
    i=i+1
    print(i)
}

USED_CUTOFF=1

summary(NUM)
tiff('GRID3.tiff',width=9,height=9, res=400,units='in')
COLS=rep('grey80',length(NUM))#colorRampPalette(c("grey80","yellow",'red1',"red3"))(max(NUM))
COLS[which(NUM>=USED_CUTOFF)]='yellow3'
COLS[which(NUM>=USED_CUTOFF+1)]='red3'
plot(V_CELL_VEC,pch=15,col=COLS,cex=0.4)
dev.off()


tiff('GRID4.tiff',width=9,height=9, res=400,units='in')
plot(V_CELL_VEC[which(NUM>=USED_CUTOFF),],pch=15,col='grey50',cex=0.4)
dev.off()

USE_VEC=V_CELL_VEC[which(NUM>=USED_CUTOFF),]
#plot(density(NUM))
saveRDS(USE_VEC,file='USE_VEC.RDS')
saveRDS(CELL_LABEL,file='CELL_LABEL.RDS')
saveRDS(V_CELL_VEC, file='V_CELL_VEC.RDS')
saveRDS(NUM,file='NUM.RDS')

#########################################
##########################################


V_CELL_VEC=readRDS(file='V_CELL_VEC.RDS')
CELL_LABEL=readRDS(file='CELL_LABEL.RDS')
META=readRDS('META.RDS')

BATCH=names(table(META$orig.ident))

M=c()
i=1
while(i<=nrow(V_CELL_VEC)){
    #this_vec=V_CELL_VEC[i,]
    this_name=rownames(V_CELL_VEC)[i]
    this_index=which(CELL_LABEL ==this_name)
    this_m=table(META$orig.ident[this_index])
    
    M=cbind(M,this_m)   
    i=i+1
    print(i)
}

M=t(M)
rownames(M)=rownames(V_CELL_VEC)

saveRDS(M,file='M.RDS')
############################################










