G3=readRDS('G3.RDS')
G4=readRDS('G4.RDS')

dim(G3)
dim(G4)

DATA=cbind(G3,G4)

rm(G3)
rm(G4)
gc()


library(Seurat)

G1=as.character(read.table('./MARKER/G3MYC',header=FALSE)[,1])
G2=as.character(read.table('./MARKER/G3NRL',header=FALSE)[,1])
G3=as.character(read.table('./MARKER/G4UD',header=FALSE)[,1])
G4=as.character(read.table('./MARKER/G4',header=FALSE)[,1])


MKDATA=DATA[which(rownames(DATA) %in% c(G1,G2,G3,G4)),]
pbmc <- CreateSeuratObject(counts = MKDATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


EXP1=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G1),]
EXP2=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G2),]
EXP3=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G3),]
EXP4=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G4),]


S1=apply(EXP1,2,mean)
S2=apply(EXP2,2,mean)
S3=apply(EXP3,2,mean)
S4=apply(EXP4,2,mean)


SS1=pnorm(scale(S1))
SS2=pnorm(scale(S2))
SS3=pnorm(scale(S3))
SS4=pnorm(scale(S4))

P1=c(-1,1)
P2=c(1,1)
P3=c(-1,-1)
P4=c(1,-1)


X=(SS1*P1[1]+SS2*P2[1]+SS3*P3[1]+SS4*P4[1])/4
Y=(SS1*P1[2]+SS2*P2[2]+SS3*P3[2]+SS4*P4[2])/4


############
VariableFeatures(object = pbmc)=all.genes
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)



G3=c('BT2016012G3','BT2018034G3','BT309','BT311','BT316','BT317','BT319')
G4=c('BT2016040G4','BT2016137G4','BT2017035G4','BT2017075G4','BT308','MB3076','MB3271','BT318')
SHH=c('BT2017017SHH','SM4217SHH','BT312')
EPN=c('BT2016062PFA','BT2016090PFA8','BT2017081PFA','BT2018022PFA9','BT307','BT314','BT315')
ATRT=c('BT313')
CPA=c('BT2017065PA','PA0706','PA2406')


pbmc@meta.data$type=rep('NA',ncol(pbmc))
pbmc@meta.data$type[which(Idents(pbmc) %in% G3)]='G3'
pbmc@meta.data$type[which(Idents(pbmc) %in% G4)]='G4'



pbmc@reductions$umap@cell.embeddings[,1]=X
pbmc@reductions$umap@cell.embeddings[,2]=Y

######
saveRDS(pbmc,'pbmc_tiny_G3G4.RDS')
######


##################################################################
##################################################################
##################################################################
#Local
library('Seurat')
pbmc=readRDS('pbmc_tiny_G3G4.RDS')

DimPlot(pbmc)



G1=as.character(read.table('./MARKER/G3MYC',header=FALSE)[,1])
G2=as.character(read.table('./MARKER/G3NRL',header=FALSE)[,1])
G3=as.character(read.table('./MARKER/G4UD',header=FALSE)[,1])
G4=as.character(read.table('./MARKER/G4',header=FALSE)[,1])


#MKDATA=DATA[which(rownames(DATA) %in% c(G1,G2,G3,G4)),]
pbmc <- CreateSeuratObject(counts = MKDATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#source('BEER.R')
#EXP=as.matrix(pbmc@assays$RNA@data)
#BATCH=as.character(pbmc@meta.data$orig.ident)
#NEWEXP=.combat(EXP,BATCH)
#pbmc@assays$RNA@data=NEWEXP

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


EXP1=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G1),]
EXP2=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G2),]
EXP3=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G3),]
EXP4=pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% G4),]

#NEXP1=t(apply(EXP1,1,pnorm))
#NEXP2=t(apply(EXP2,1,pnorm))
#NEXP3=t(apply(EXP3,1,pnorm))
#NEXP4=t(apply(EXP4,1,pnorm))



S1=apply(EXP1,2,mean)
S2=apply(EXP2,2,mean)
S3=apply(EXP3,2,mean)
S4=apply(EXP4,2,mean)

S1=scale(S1)
S2=scale(S2)
S3=scale(S3)
S4=scale(S4)

S1=S1-min(S1)
S2=S2-min(S2)
S3=S3-min(S3)
S4=S4-min(S4)

S1=S1/max(S1)
S2=S2/max(S2)
S3=S3/max(S3)
S4=S4/max(S4)

#S1=pnorm(S1)
#S2=pnorm(S2)
#S3=pnorm(S3)
#S4=pnorm(S4)


BIND=cbind(S1,S2,S3,S4)
VAR=apply(BIND,1,var)

SS1=S1
SS2=S2
SS3=S3
SS4=S4


SS1=S1*log(VAR+1,10)
SS2=S2*log(VAR+1,10)
SS3=S3*log(VAR+1,10)
SS4=S4*log(VAR+1,10)


P1=c(-1,1)
P2=c(1,1)
P3=c(-1,-1)
P4=c(1,-1)


X=(SS1*P1[1]+SS2*P2[1]+SS3*P3[1]+SS4*P4[1])
Y=(SS1*P1[2]+SS2*P2[2]+SS3*P3[2]+SS4*P4[2])
X=scale(X)
Y=scale(Y)



pbmc@reductions$umap@cell.embeddings[,1]=X
pbmc@reductions$umap@cell.embeddings[,2]=Y

DimPlot(pbmc)

DimPlot(pbmc,group.by='type')
FeaturePlot(pbmc,features=c('MYC','NRL','SOX11','ERF'))




