
setwd('C:/Users/cchmc/Desktop/BEER')
library(liger)
library(cowplot)



D1=readRDS('MGH36.RDS')
D2=readRDS('MGH53.RDS')
D3=readRDS('MGH54.RDS')
D4=readRDS('MGH60.RDS')
D5=readRDS('MGH93.RDS')
D6=readRDS('MGH97.RDS')

pbmc.data = list(D1=D1,D2=D2,D3=D3,D4=D4,D5=D5,D6=D6)

# Create liger object
a.pbmc <- createLiger(pbmc.data)
a.pbmc <- normalize(a.pbmc)
# Can pass different var.thresh values to each dataset if one seems to be contributing significantly
# more genes than the other
a.pbmc <- selectGenes(a.pbmc, var.thresh = c(0.3, 0.875), do.plot = F)
#s.var.genes <- #readRDS('~/Downloads/pbmc_alignment/var_genes.RDS')
#a.pbmc@var.genes <- s.var.genes
a.pbmc <- scaleNotCenter(a.pbmc)
#k.suggest <- suggestK(a.pbmc, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F,
#                     nrep = 5)
a.pbmc <- optimizeALS(a.pbmc, k=22, thresh = 5e-5, nrep = 3)

a.pbmc <- runTSNE(a.pbmc, use.raw = T)
p1 <- plotByDatasetAndCluster(a.pbmc, return.plots = T)
# Plot by dataset
#print(p1[[1]])

a.pbmc <- quantileAlignSNF(a.pbmc, resolution = 0.4, small.clust.thresh = 20)

a.pbmc <- runTSNE(a.pbmc)
p_a <- plotByDatasetAndCluster(a.pbmc, return.plots = T) 
# Modify plot output slightly
p_a[[1]] <- p_a[[1]] + theme_classic() + theme(legend.position = c(0.85, 0.15)) + 
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))

print(p_a[[1]])






###############################################################################3


setwd('C:/Users/cchmc/Desktop/BEER')
library(liger)
library(cowplot)


SEED=135
T.NUM=2
C.NUM=300
G.NUM=2000
P.MEAN=0
R.SIZE=0.1
A.MEAN=50
#DR=0.1 # Drop-out


set.seed(SEED)
DATA=matrix(0,ncol=T.NUM*C.NUM,nrow=G.NUM)
rownames(DATA)=paste0('G',1:nrow(DATA))
TMP1=rep(paste0('T',1:T.NUM),time=1,each=C.NUM)
TMP2=rep(1:C.NUM,time=T.NUM)
colnames(DATA)=paste0(TMP1,'_',TMP2)


this_start=rpois(G.NUM, P.MEAN)
j=0
while(j<T.NUM){
    
    this_a=rnorm(G.NUM)
    #this_b=rnorm(G.NUM)
  
    i=1
    while(i<= C.NUM){
        this_col=C.NUM*j+i
        this_r=rnorm(G.NUM)*R.SIZE
        this_add=rpois(G.NUM,A.MEAN)
        
        this_data=this_start+(this_a+this_r)*i+this_add
        
        this_data[which(this_data<0)]=0
        this_data=round(this_data)
        #this_data[sample(1:G.NUM,round(DR*G.NUM))]=0      
        DATA[,this_col]=this_data
        i=i+1}

    j=j+1}


saveRDS(DATA,'RDATA.RDS')
########################################


library(Seurat)

pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
VariableFeatures(object = pbmc)=all.genes
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
#pbmc <- RunUMAP(pbmc, dims = 1:10)
#DimPlot(pbmc, reduction = "umap")
pbmc <- RunUMAP(pbmc, dims = 1:10,n.components = 2)
umap=pbmc@reductions$umap@cell.embeddings

require(scales)
my_color_palette <- hue_pal()(T.NUM)
COL=rep(my_color_palette,time=1,each=C.NUM)
plot(umap[,1],umap[,2],col=COL,pch=16,)





ligerex = createLiger(list(D1=DATA)) #Can also pass in more than 2 datasets
ligerex = normalize(ligerex)
ligerex = selectGenes(ligerex, var.thresh = 0.1)
ligerex = scaleNotCenter(ligerex)
ligerex = optimizeALS(ligerex, k = 10) 
ligerex = quantileAlignSNF(ligerex)
ligerex = runUMAP(ligerex,k=2)
#plotByDatasetAndCluster(ligerex)

umap=ligerex@tsne.coords
require(scales)
my_color_palette <- hue_pal()(T.NUM)
COL=rep(my_color_palette,time=1,each=C.NUM)
plot(umap[,1],umap[,2],col=COL,pch=16,)



























#############
#3D
ligerex = runUMAP(ligerex,k=3)
library("rgl")
library("car")
scatter3d(umap[,1], umap[,2], umap[,3], point.col = COL, surface=FALSE)
umap=ligerex@tsne.coords
require(scales)
my_color_palette <- hue_pal()(T.NUM)
COL=rep(my_color_palette,time=1,each=C.NUM)
library("rgl")
library("car")
scatter3d(umap[,1], umap[,2], umap[,3], point.col = COL, surface=FALSE)

saveRDS(ligerex,'LIGER.RDS')

####################################################





pbmc <- RunUMAP(pbmc, dims = 1:10,n.components = 3)
umap=pbmc@reductions$umap@cell.embeddings

require(scales)

my_color_palette <- hue_pal()(T.NUM)
COL=rep(my_color_palette,time=1,each=C.NUM)
library("rgl")
library("car")
scatter3d(umap[,1], umap[,2], umap[,3], point.col = COL, surface=FALSE)





