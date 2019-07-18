setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')


REF=.readTable(PATH='GSE92332_atlas_UMIcounts.txt',SEP='\t')


getbatch <- function(x){
    y=unlist(strsplit(x, "_"))
    y=y[length(y)]
    return(y)
}
CN=colnames(REF)
BATCH=apply(matrix(CN,ncol=1),1,getbatch)
table(BATCH)
LABEL=BATCH

library(dplyr)
library(Seurat)

#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO<- Read10X(data.dir = "./Small_Intestine_KO")
AGE <- Read10X(data.dir = "./age")
YOUNG <- Read10X(data.dir = "./young")

BATCH=c(rep('NATURE',ncol(REF)),
        rep('CDC42HET',ncol(CDC42HET)),
        rep('CDC42KO',ncol(CDC42KO)),   
        rep('AGE',ncol(AGE)),
       rep('YOUNG',ncol(YOUNG))
       )
####################

D1=.simple_combine(REF, CDC42HET)$combine
rm(REF)
rm(CDC42HET)
gc()


D2=.simple_combine(CDC42KO, AGE)$combine
rm(CDC42KO)
rm(AGE)
gc()


D3=.simple_combine(D1, D2)$combine
rm(D1)
rm(D2)
gc()

DATA=.simple_combine(D3, YOUNG)$combine
rm(D3)
rm(YOUNG)
gc()
saveRDS(DATA,'DATA.RDS')
saveRDS(BATCH,'BATCH.RDS')


TMP=rep('NA',ncol(DATA))
TMP[which(BATCH=='NATURE')]=LABEL
LABEL=TMP
saveRDS(BATCH,'LABEL.RDS')


###################################################################

getNon0=function(x){
   return(length(which(x>0)))
   }

NON0=apply(DATA,2,getNon0)

USED.CELL=which(NON0>800 & NON0<6000)


DATA.USED=DATA[,USED.CELL]
BATCH.USED=BATCH[USED.CELL]
LABEL.USED=LABEL[USED.CELL]
rm(DATA)
gc()


NON0.G=apply(DATA.USED,1,getNon0)
USED.G=which(NON0.G>50)
DATA.USED=DATA.USED[USED.G,]
dim(DATA.USED)
#[1] 13913 28866


saveRDS(DATA.USED,'DATA.USED.RDS')
saveRDS(BATCH.USED,'BATCH.USED.RDS')
saveRDS(LABEL.USED,'LABEL.USED.RDS')

###################################################################


mybeer=BEER(DATA.USED, BATCH.USED, GNUM=30, PCNUM=50, ROUND=1, CPU=2, GN=3000, SEED=1, MTTAG='^mt-',REGBATCH=TRUE)   
saveRDS(mybeer,'mybeer.RDS')

####################################################

# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


pbmc <- mybeer$seurat
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)    


pbmc <- mybeer$seurat
PCUSE <- mybeer$select
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)

DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)   


pbmc <- mybeer$seurat
PCUSE=mybeer$select
pbmc=BEER.combat(pbmc) 
umap=BEER.bbknn(pbmc, PCUSE, NB=4, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)




####################################################



pbmc@meta.data$celltype=LABEL.USED
pbmc@meta.data$celltype[LABEL.USED=='NA']=NA

DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=F)

DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)

####################################################





#######
# Transfer
VEC=pbmc@reductions$umap@cell.embeddings
set.seed(123)
N=500
K=kmeans(VEC,centers=N)
pbmc@meta.data$kclust=K$cluster   
#DimPlot(pbmc, reduction.use='umap', group.by='kclust', pt.size=0.1,label=T)

pbmc@meta.data$transfer=rep(NA, length(pbmc@meta.data$celltype))
TMP=cbind(pbmc@meta.data$celltype, pbmc@meta.data$kclust)

KC=unique(pbmc@meta.data$kclust)
i=1
while(i<=length(KC)){
    this_kc=KC[i]
    this_index=which(pbmc@meta.data$kclust==this_kc)
    this_tb=table(pbmc@meta.data$celltype[this_index])
    if(length(this_tb)!=0){
        this_ct=names(this_tb)[which(this_tb==max(this_tb))[1]]
        pbmc@meta.data$transfer[this_index]=this_ct}
    i=i+1}
    
pbmc@meta.data$tf.ct=pbmc@meta.data$celltype
NA.index=which(is.na(pbmc@meta.data$celltype))
pbmc@meta.data$tf.ct[NA.index]=pbmc@meta.data$transfer[NA.index]



RNA.cells=colnames(pbmc)[which(pbmc@meta.data$batch=='NATURE')]
ATAC.cells=colnames(pbmc)[which(pbmc@meta.data$batch!='NATURE')]




library(ggplot2)

plot.all <- DimPlot(pbmc, reduction.use='umap', group.by='batch', 
    pt.size=0.1,label=F) + labs(title = "Batches")

plot.ct <- DimPlot(pbmc,reduction.use='umap', group.by='tf.ct', 
    pt.size=0.1,label=T) + labs(title = "CellType")

plot.rna <- DimPlot(pbmc, cells=RNA.cells,reduction.use='umap', 
    group.by='tf.ct', pt.size=0.1,label=T) + labs(title = "Nature")

plot.atac <- DimPlot(pbmc, cells=ATAC.cells,reduction.use='umap', 
    group.by='tf.ct', pt.size=0.1,label=T) + labs(title = "OurData")

CombinePlots(list(all=plot.all, ct=plot.ct, rna=plot.rna, atac=plot.atac))


saveRDS(pbmc,'pbmc_final.RDS')

TAB=table(pbmc@meta.data$tf.ct, pbmc@meta.data$batch)
.writeTable(DATA=TAB,PATH='TABLE.txt')

SUM=apply(TAB,2,sum)

TAB.100=round(TAB/SUM*100)
.writeTable(DATA=TAB.100,PATH='TABLE100.txt')


FeaturePlot(pbmc, features = c('Kit','Cd4'))
