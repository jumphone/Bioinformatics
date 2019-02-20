
library(Seurat)
library(dplyr)
source('scRef.R')

####################################
pbmc=readRDS('13Nature_10X.RDS')
MAP=read.table('Clusters_combine.txt',sep='\t',header=T)

MAP[,1]=as.character(MAP[,1])
MAP[,2]=as.character(MAP[,2])
TAG=pbmc@meta.data$C
TAG=as.character(TAG)

tmp=c()
i=1
while(i<=length(TAG)){
   if(TAG[i] %in% MAP[,2]){
       this_tmp=MAP[which(MAP[,2]==TAG[i]),1]
   }else{this_tmp='NA'}
   tmp=c(tmp,this_tmp)
   i=i+1
}

ref_tag = cbind(colnames(pbmc@data),tmp)

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]

NatureRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )
saveRDS(NatureRef,file='NatureRef.RDS')

##################################################

pbmc.data.1=read.table('decidua-20181218_exon_tagged.dge.txt',sep='\t',header=T,row.names=1)
pbmc.data.2=read.table('placenta-20181218_exon_tagged.dge.txt',sep='\t',header=T,row.names=1)
pbmc.1 <- CreateSeuratObject(raw.data = pbmc.data.1, min.cells = 3, min.genes = 200, project = "decidua")
pbmc.2 <- CreateSeuratObject(raw.data = pbmc.data.2, min.cells = 3, min.genes = 200, project = "placenta")

pbmc=MergeSeurat(pbmc.1, pbmc.2, add.cell.id1='decidua', add.cell.id2='placenta')


mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.4))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F, mean.function = ExpMean, dispersion.function = LogVMR,  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito","orig.ident"))

PCNUM=200
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, pcs.compute=PCNUM, genes.print = 5)
PCElbowPlot(object = pbmc,num.pc=PCNUM)

PCUSE=1:24
#pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)
pbmc <- RunUMAP(pbmc, dims.use = PCUSE)

#TSNEPlot(object = pbmc)
DimPlot(pbmc, reduction.use = "umap", pt.size=0.2)

TSNE_VEC=pbmc@dr$umap@cell.embeddings
D=dist(TSNE_VEC)
H=hclust(D)
C=cutree(H, k=12) 
pbmc@meta.data$C=C


DimPlot(pbmc, reduction.use = "umap", pt.size=2, group.by='orig.ident')
DimPlot(pbmc, reduction.use = "umap", pt.size=2, group.by='C')



##################################################

source('scRef.R')

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

LocalAgg=.generate_ref(exp_sc_mat, cbind(names(pbmc@ident),as.character(pbmc@meta.data$C)), min_cell = 1 )

out=.get_cor(LocalAgg, NatureRef, method='kendall',CPU=4, print_step=10)
tag=.get_tag_max(out)

tag_all=c()
for(one in pbmc@meta.data$C){
   tag_all=c(tag_all, tag[which(tag[,1]==as.character(one)),2] )
}

#out=SCREF(exp_sc_mat,NatureRef,min_cell=10)
#tag_all=out$tag2[,2]

pbmc@meta.data$scref=as.character(tag_all)
DimPlot(pbmc, reduction.use = "umap",group.by='scref', do.label=T, pt.size=2, label.size=5)

saveRDS(pbmc,file='HongFangZi.RDS')

pdf('HongFangZi.pdf',width=7,height=5)
DimPlot(pbmc, reduction.use = "umap", pt.size=2, do.label=T, group.by='orig.ident')
DimPlot(pbmc, reduction.use = "umap", pt.size=2, do.label=T, group.by='C')
DimPlot(pbmc, reduction.use = "umap",group.by='scref', do.label=T, pt.size=2, label.size=5)
dev.off()

tmp=pbmc@ident
pbmc@ident=as.factor(pbmc@meta.data$scref)
names(pbmc@ident)=names(tmp)


pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
    thresh.use = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)


pdf('HongFangZi_withmarker.pdf',width=15,height=12)
DimPlot(pbmc, reduction.use = "umap", pt.size=2, do.label=T, group.by='orig.ident')
DimPlot(pbmc, reduction.use = "umap", pt.size=2, do.label=T, group.by='C')
DimPlot(pbmc, reduction.use = "umap",group.by='scref', do.label=T, pt.size=2, label.size=5)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()


saveRDS(pbmc.markers,'HongFangZi_marker.RDS')

G102=read.table('102G.txt',sep='\t')

exp_data=as.matrix(pbmc@data)

#exp_data=as.matrix(pbmc@scale.data)

#LIM=1.5
#exp_data[which(exp_data>LIM)]=LIM
#exp_data[which(exp_data< -LIM)]= -LIM

V=which(rownames(exp_data) %in% as.character(G102[,1]))
O=order(pbmc@meta.data$scref)
OT=pbmc@meta.data$scref[O]
data=as.matrix(exp_data[V,O])
CLUSTER_NUM=length(unique(pbmc@meta.data$scref))


CCC=rainbow(CLUSTER_NUM)
COL=rep('red',length(data[1,]))
tmp=''
i=0
j=1
for(one in OT){
   if(one !=tmp){i=i+1}
   COL[j]=CCC[i]
   tmp=one
   j=j+1
}


pdf('HEATMAP.pdf',width=10,height=12)
heatmap.2(data,scale=c("none"),labCol = F,dendrogram='none',Colv=F,trace='none',ColSideColors=COL,
          col=colorRampPalette(c('blue3','grey95','red3')),main=as.character(i) ,margins=c(10,10))

plot(1:length(unique(OT)),col=CCC,pch=16,cex=3)
text(1:length(unique(OT)),labels=unique(OT),pos=3)
dev.off()



######################################################

#install.packages("ggcorrplot")
library('ggcorrplot')
tmp=cor(as.matrix(pbmc@data))
ggcorrplot(tmp[1:10,1:10], method = "circle")

LR=read.table('LR.txt',sep='\t')
GENE=rownames(pbmc@data)
V=c()
i=1
while(i<=length(LR[,1])){
if(LR[i,1] %in% GENE && LR[i,2] %in% GENE){V=c(V,i)}
i=i+1
}

VLR=LR[V,]

library('survcomp')
LRG=c(as.character(VLR[,1]),as.character(VLR[,2]))

TAG=pbmc@meta.data$scref
EXP=pbmc@data
G=rownames(pbmc@data)

EXPV=as.matrix(EXP[which(G %in% LRG),])

UT=unique(TAG)


PVM=c()
for(one in UT){
 V=which(TAG==one)
 N=which(TAG!=one)
 pvs=c()
 i=1
 while(i <=nrow(EXPV)){
   #gene=rownames(EXPV)[i]
   pv=wilcox.test(EXPV[i,V] , EXPV[i,N], alternative='greater')$p.value
   pvs=c(pvs,pv)
   i=i+1}
 PVM=cbind(PVM,pvs)
}

colnames(PVM)=UT
rownames(PVM)=rownames(EXPV)

VLR[,1]=as.character(VLR[,1])
VLR[,2]=as.character(VLR[,2])


OUT=c()
i=1
while(i<=nrow(VLR)){
   
    R=VLR[i,1]
    L=VLR[i,2]
   
    dNKR=PVM[which(rownames(PVM)==R),which(colnames(PVM)=='dNK')]
    dNKL=PVM[which(rownames(PVM)==L),which(colnames(PVM)=='dNK')]
   
    tmp=c()
    for(one in UT){
        thisL=PVM[which(rownames(PVM)==L),which(colnames(PVM)==one)]
        thisR=PVM[which(rownames(PVM)==R),which(colnames(PVM)==one)]
        dNKR_P=combine.test(c(dNKR, thisL), method ='fisher')
        dNKL_P=combine.test(c(dNKL, thisR), method ='fisher')
        tmp=c(tmp,dNKR_P,dNKL_P)
    }
   
   OUT=cbind(OUT,tmp)
   
    i=i+1}

tmp1=paste0(UT,'_2')
tmp2=paste0(UT,'_1')



rownames(OUT)=as.character(t(cbind(tmp1,tmp2)))



colnames(OUT)=paste0(VLR[,1],'_',VLR[,2])

#NEWOUT=apply(OUT,2,p.adjust,'fdr')
NEWOUT=OUT
NEWOUT=-log(NEWOUT,10)
NEWOUT=apply(NEWOUT,2,scale)

NEWOUT[which(NEWOUT>2)]=2
NEWOUT[which(NEWOUT< -2)]=-2

colnames(NEWOUT)=colnames(OUT)
rownames(NEWOUT)=rownames(OUT)

library('gplots')
heatmap.2(t(NEWOUT),scale=c("none"),dendrogram='none',Colv=F,trace='none',col=colorRampPalette(c('royalblue','grey95','indianred')) ,margins=c(10,10))





