
#setwd('F:/Zhenglab/Combine')

setwd('/Volumes/Feng/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#REF=.readTable(PATH='GSE92332_atlas_UMIcounts.txt',SEP='\t')
REF=readRDS('REF.RDS')
getbatch <- function(x){
    y=unlist(strsplit(x, "_"))
    y=y[length(y)]
    return(y)
}

#saveRDS(REF,file='REF.RDS')


CN=colnames(REF)
BATCH=apply(matrix(CN,ncol=1),1,getbatch)
table(BATCH)
LABEL=BATCH



############################
REF=readRDS('REF.RDS')
getbatch <- function(x){
    y=unlist(strsplit(x, "_"))
    y=y[length(y)]
    return(y)
}

#saveRDS(REF,file='REF.RDS')
CN=colnames(REF)
BATCH=apply(matrix(CN,ncol=1),1,getbatch)
table(BATCH)
LABEL=BATCH
saveRDS(LABEL,file='LABEL.RDS')

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
LocalRef=.generate_ref(REF, cbind(LABEL,LABEL),M='MEAN', min_cell=1)  

saveRDS(LocalRef,file='LocalRef.RDS')

############################

CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO<- Read10X(data.dir = "./Small_Intestine_KO")
AGE <- Read10X(data.dir = "./age")
YOUNG <- Read10X(data.dir = "./young")

BATCH1=c(rep('NATURE',ncol(REF)),
        rep('CDC42HET',ncol(CDC42HET)),
        rep('CDC42KO',ncol(CDC42KO))
       )

BATCH2=c(rep('NATURE',ncol(REF)),
        rep('AGE',ncol(AGE)),
       rep('YOUNG',ncol(YOUNG))
)
####################

D1=.simple_combine(REF, CDC42HET)$combine
rm(CDC42HET)
gc()


DATA1=.simple_combine(D1, CDC42KO)$combine
rm(CDC42KO)
gc()

D2=.simple_combine(REF, AGE)$combine
rm(AGE)
gc()


DATA2=.simple_combine(D2, YOUNG)$combine
rm(YOUNG)
gc()

####################
saveRDS(DATA1,'DATA1.RDS') # Zheng Zhang
saveRDS(BATCH1,'BATCH1.RDS')

saveRDS(DATA1,'DATA2.RDS') # Age Young
saveRDS(BATCH1,'BATCH2.RDS')
####################


#######################

setwd('/Volumes/Feng/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
.set_python('/Users/zha8dh/anaconda3/bin/python')

DATA=readRDS('DATA1.RDS') # Zheng Zhang
BATCH=readRDS('BATCH1.RDS')


mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=3, GN=5000, SEED=1, COMBAT=TRUE)
saveRDS(mybeer,file='mybeer1.RDS')




PCUSE=mybeer$select
#PCUSE=.selectUSE(mybeer, CUTR=0.8, CUTL=0.8, RR=0.5, RL=0.5)

COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

pbmc <- mybeer$seurat  
PCUSE=mybeer$select
#PCUSE=.selectUSE(mybeer, CUTR=0.8, CUTL=0.8, RR=0.5, RL=0.5)

#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)

saveRDS(pbmc,file='pbmc1.RDS')

FeaturePlot(pbmc,features=c('Ptprc'))





########

LABEL=readRDS(file='LABEL.RDS')
pbmc@meta.data$celltype=rep(NA,ncol(pbmc))
pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='NATURE')]=LABEL
#DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=F)
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)





############################
############################
############################






setwd('/Volumes/Feng/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
.set_python('/Users/zha8dh/anaconda3/bin/python')

DATA=readRDS('DATA1.RDS') # Zheng Zhang
BATCH=readRDS('BATCH1.RDS')


pbmc=readRDS('./pbmc1.RDS')


LABEL=readRDS(file='LABEL.RDS')
pbmc@meta.data$celltype=rep(NA,ncol(pbmc))
pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='NATURE')]=LABEL
#DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=F)
DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)




#Stem cell: "Lgr5, Ascl2, Gkn3, Slc12a2, Axin2, Olfm4". 
#For the TA (cycling) cells, please check "Mki67, Cdk4, Mcm5, Mcm6, Pcna", 
#Enterrocyte markers are "Alpi, Apoa1, Apoa4, Fabp1", 
#Goblet cell markers are "Muc2, Clca3, Tff3, Agr2". 
#Paneth cell markers are" Lyz1, Defa17, Defa22, Defa24, Ang4". 
#Endocrine cell markers are "Chga, Chgb, Tac1, Tph1, Neurog3"

USED_CELL=which(pbmc@meta.data$batch %in% c('CDC42KO','CDC42HET'))

pbmc_zhengzhang <- CreateSeuratObject(counts = pbmc@assays$RNA@counts[,USED_CELL], project = "ZhengZhang", min.cells = 0, min.features = 0)
pbmc_zhengzhang@meta.data=pbmc@meta.data[USED_CELL,]
pbmc_zhengzhang <- NormalizeData(pbmc_zhengzhang, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc_zhengzhang)
pbmc_zhengzhang <- ScaleData(pbmc_zhengzhang, features = all.genes)
VariableFeatures(pbmc_zhengzhang)=VariableFeatures(pbmc)
pbmc_zhengzhang <- RunPCA(pbmc_zhengzhang, npcs=10,features = VariableFeatures(object = pbmc_zhengzhang))
pbmc_zhengzhang <- RunUMAP(pbmc_zhengzhang, dims = 1:10)

pbmc_zhengzhang@reductions$umap@cell.embeddings=pbmc@reductions$umap@cell.embeddings[USED_CELL,]



pdf('/Users/zha8dh/Desktop/CCHMC_Project/ZhengLab/FIG/ZhengZhangCells.pdf')
this_pbmc=pbmc_zhengzhang
FeaturePlot(this_pbmc,  features=c('Lgr5', 'Ascl2', 'Gkn3', 'Slc12a2', 'Axin2', 'Olfm4')) #Stem cell
FeaturePlot(this_pbmc,  features=c('Mki67', 'Cdk4', 'Mcm5', 'Slc12a2', 'Mcm6', 'Pcna')) #TA (cycling) cells
FeaturePlot(this_pbmc,  features=c('Alpi', 'Apoa1', 'Apoa4', 'Fabp1')) #Enterrocyte
FeaturePlot(this_pbmc,  features=c('Muc2', 'Clca3', 'Tff3', 'Agr2')) # Goblet
FeaturePlot(this_pbmc,  features=c('Lyz1', 'Defa17', 'Defa22', 'Defa24','Ang4'))#Paneth
FeaturePlot(this_pbmc,  features=c('Chga', 'Chgb', 'Tac1', 'Tph1', 'Neurog3'))#Endocrine
dev.off()

pdf('/Users/zha8dh/Desktop/CCHMC_Project/ZhengLab/FIG/AllCells.pdf')
this_pbmc=pbmc
FeaturePlot(this_pbmc,  features=c('Lgr5', 'Ascl2', 'Gkn3', 'Slc12a2', 'Axin2', 'Olfm4')) #Stem cell
FeaturePlot(this_pbmc,  features=c('Mki67', 'Cdk4', 'Mcm5', 'Slc12a2', 'Mcm6', 'Pcna')) #TA (cycling) cells
FeaturePlot(this_pbmc,  features=c('Alpi', 'Apoa1', 'Apoa4', 'Fabp1')) #Enterrocyte
FeaturePlot(this_pbmc,  features=c('Muc2', 'Clca3', 'Tff3', 'Agr2')) # Goblet
FeaturePlot(this_pbmc,  features=c('Lyz1', 'Defa17', 'Defa22', 'Defa24','Ang4'))#Paneth
FeaturePlot(this_pbmc,  features=c('Chga', 'Chgb', 'Tac1', 'Tph1', 'Neurog3'))#Endocrine
dev.off()




############################
############################
############################






#######
# Transfer
VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)
N=100
K=kmeans(VEC,centers=N)


pbmc@meta.data$kclust=K$cluster   

FeaturePlot(pbmc,features=c('Ptprc'))


DimPlot(pbmc, reduction.use='umap', group.by='kclust', pt.size=0.1,label=T) + NoLegend()


#pbmc@meta.data$kclust[which(pbmc@meta.data$kclust %in% c(55,152,53,33,84,30,93,36,157,118,120,47,126,57))]='Imm'
pbmc@meta.data$kclust[which(pbmc@meta.data$kclust %in% c(89,79,66,18,76,4,54,45,99,32,25))]='Imm'

DimPlot(pbmc, reduction.use='umap', group.by='kclust', pt.size=0.1,label=T) + NoLegend()

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
LocalCLUST=.generate_ref(as.matrix(pbmc@assays$RNA@data), cbind(pbmc@meta.data$kclust,pbmc@meta.data$kclust),M='MEAN', min_cell=1)  

LocalCLUST=LocalCLUST[,which(colnames(LocalCLUST)!='Imm')]


OUT=SCREF(LocalCLUST, LocalRef, CPU=4, print_step=10,min_cell=1)
TAG=OUT$tag2



pbmc@meta.data$scref=pbmc@meta.data$kclust
i=1
while(i<=nrow(TAG)){
    this_index=which(as.character(pbmc@meta.data$kclust)==TAG[i,1])
    pbmc@meta.data$scref[this_index]=TAG[i,2]
    i=i+1}





DimPlot(pbmc, reduction.use='umap', group.by='scref', pt.size=0.1,label=T)



saveRDS(pbmc,file='pbmc1.RDS')









##################
mybeer2=BEER(DATA2, BATCH2, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE)












