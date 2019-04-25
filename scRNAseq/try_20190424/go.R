library(Seurat)
library(dplyr)
library(Matrix)
#pbmc.data <- read.table("MGH54_mat.txt", sep='\t',row.names=1, header=T)

load('Seurat_EXP_cluster.Robj')
pbmc.raw.data=as.matrix(EXP_cluster@raw.data[,which(colnames(EXP_cluster@raw.data) %in% colnames(EXP_cluster@scale.data))])
pbmc.data=as.matrix(EXP_cluster@scale.data)
#pbmc.data=as.matrix(EXP_cluster@data)

#LR=read.table('RL.txt',header=T,sep='\t')
LR=read.table('RL_mouse.txt',header=T,sep='\t')

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

rm(EXP_cluster)
gc()

ONE=.data2one(pbmc.raw.data, rownames(pbmc.data), CPU=4, PCNUM=50, SEED=123,  PP=30)

WINDOW=300

RANK=rank(ONE)
LENGTH=length(ONE)
BIN=c()
i=1
while(i<=LENGTH/WINDOW  ){
this_index=which((i-1)*WINDOW< RANK & i*WINDOW>=RANK)
BIN=cbind(BIN,this_index)
i=i+1
}

EXP=pbmc.data
GENE=rownames(EXP)

saveRDS(BIN,'BIN.RDS')
saveRDS(ONE,'ONE.RDS')
################################################

ALL=c(as.character(LR[,1]),as.character(LR[,2]) )

permu_gene_index=which(GENE %in% ALL)

set.seed(123)
TIME=10000
MEAN=EXP[permu_gene_index,c(1:TIME)]*0
colnames(MEAN)=as.character(c(1:TIME))
i=1
while(i<=TIME){
this_index=sample(ncol(EXP),WINDOW)
this_mean=apply(EXP[permu_gene_index,this_index],1,mean)
MEAN[,i]=this_mean
i=i+1
if(i%%100==1){print(i)}}

saveRDS(MEAN,file=paste0('MEAN.RDS' ))

################################################
EXP_LR=EXP[permu_gene_index,]


ECDF=c()
j=1
while(j<=nrow(MEAN)){
ECDF=c(ECDF,ecdf(as.numeric(MEAN[j,])))
j=j+1}


PMAT = MEAN[,c(1:ncol(BIN))]*0
i=1
while(i<=ncol(BIN)){
this_bin_index=BIN[,i]
  
this_bin_mean_exp=apply(EXP_LR[,this_bin_index],1,mean)  
this_p_list=c()
j=1
while(j<=length(this_bin_mean_exp)){
  this_p=-log(1+1/TIME-ECDF[[j]](this_bin_mean_exp[j]),10)
  this_p_list=c(this_p_list,this_p)
  j=j+1
}
PMAT[,i]=this_p_list
print(i)
i=i+1
}
saveRDS(PMAT,file=paste0('PMAT.RDS' ))

################################################
CMAT=PMAT[c(1:ncol(PMAT)),]*0
rownames(CMAT)=colnames(CMAT)
rownames(CMAT)=paste0('L_',rownames(CMAT))
colnames(CMAT)=paste0('R_',colnames(CMAT))

i=1
while(i<=nrow(LR)){

this_l=LR[i,1]
this_r=LR[i,2]
if(this_l %in% GENE & this_r %in% GENE){
    this_l_index=which(rownames(PMAT)==this_l)
    this_r_index=which(rownames(PMAT)==this_r)
    this_l_bin_index=1
    while(this_l_bin_index<=nrow(CMAT)){
       this_r_bin_index=1
       while(this_r_bin_index<=ncol(CMAT)){
           CMAT[this_l_bin_index,this_r_bin_index]=CMAT[this_l_bin_index,this_r_bin_index]+ 
               PMAT[this_l_index,this_l_bin_index] - PMAT[this_r_index,this_l_bin_index] + PMAT[this_r_index,this_r_bin_index] - PMAT[this_l_index,this_r_bin_index]
        
               #PMAT[this_l_index,this_l_bin_index] - PMAT[this_r_index,this_l_bin_index] + PMAT[this_r_index,this_r_bin_index] - PMAT[this_l_index,this_r_bin_index]
               #max(PMAT[this_l_index,this_l_bin_index] - PMAT[this_r_index,this_l_bin_index],0)+ max(PMAT[this_r_index,this_r_bin_index] - PMAT[this_l_index,this_r_bin_index],0)
               #PMAT[this_l_index,this_l_bin_index]  + PMAT[this_r_index,this_r_bin_index]
           this_r_bin_index=this_r_bin_index+1
           }      
       this_l_bin_index=this_l_bin_index+1
       } 
    }
if(i%%10==1){print(i)}
i=i+1}

CMAT=as.matrix(CMAT)

saveRDS(CMAT,file=paste0('CMAT.RDS' ))
library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))

##########
#load('Seurat_EXP_cluster.Robj')

#######################
NCMAT=as.numeric(CMAT)
SNCMAT=sort(NCMAT,decreasing=T)
length(SNCMAT)*0.05

#TOP=round(length(SNCMAT)*0.05)
TOP=200
CUTOFF=SNCMAT[TOP]
#CUTOFF=SNCMAT[1]
TMP=CMAT
TMP[which(TMP<CUTOFF)]=0
PAIR=c()
SCORE=c()
i=1
while(i<=ncol(TMP)){

    if(max(TMP[,i])>0){
    this_r=i
    j=1
    while(j<=nrow(TMP)){
        if(TMP[j,i]>0){
            this_l=j
            PAIR=cbind(PAIR,c(this_l,this_r))
            SCORE=c(SCORE,TMP[j,i])}
        j=j+1}
    }

    i=i+1}

PAIR=t(PAIR)
colnames(PAIR)=c('L','R')

PAIR=PAIR[order(SCORE,decreasing=T),]
SCORE=SCORE[order(SCORE,decreasing=T)]
SCORE=round(SCORE)


#######################
#load('Seurat_EXP_cluster.Robj')
pbmc=EXP_cluster
pbmc@meta.data$bin=BIN_FLAG
VEC=pbmc@dr$tsne@cell.embeddings

BIN_FLAG=rep(NA,ncol(EXP))
i=1
while(i<=ncol(BIN)){
BIN_FLAG[BIN[,i]]=i
i=i+1
}


#library('gplots')
#heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))


par(mfrow=c(2,1))
hist(as.numeric(CMAT),xlab='SCORE',breaks=100,freq=F,main='Histogram of SCORE')
abline(v=CUTOFF,col='red')
text(x=CUTOFF,y=0,pos=4,col='red',label=as.character(round(CUTOFF)))

plot(VEC,col='grey80',pch=16,cex=0.3,main=paste0('TOP:',as.character(TOP), 
                                                 '; PERCENT:', as.character(round(TOP/length(SNCMAT)*100)),'%',
                                                 '; CUTOFF:',as.character(round(CUTOFF))))

legend("topleft", legend=c("Ligand", "Recepter"),
       fill=c("red", "blue"))

i=1
while(i<=nrow(PAIR)){
this_pair=PAIR[i,]
this_l=which(BIN_FLAG==this_pair[1])
this_r=which(BIN_FLAG==this_pair[2])
this_l_vec=VEC[this_l,]
this_r_vec=VEC[this_r,]

start_point= this_l_vec[round(runif(1)*nrow(this_l_vec)),]
end_point= this_r_vec[round(runif(1)*nrow(this_r_vec)),]
points(start_point[1],start_point[2],col='red',pch=16,cex=2)
points(end_point[1],end_point[2],col='blue',pch=16,cex=2)
points(this_l_vec,col='black',pch=16,cex=0.3)
points(this_r_vec,col='black',pch=16,cex=0.3)

segments(start_point[1], start_point[2], end_point[1],end_point[2],col='grey40',lty=3,lwd=1)
i=i+1}










##


I=1
while(I<=nrow(PAIR)){
  
png(paste0('OK/',as.character(I),'.png'),width=800,height=700)  

plot(VEC,col='grey80',pch=16,cex=0.3,main=paste0('TOP:',as.character(TOP), 
                                                 '; PERCENT:', as.character(round(TOP/length(SNCMAT)*100)),'%',
                                                 '; CUTOFF:',as.character(round(CUTOFF))))

legend("topleft", legend=c("Ligand", "Recepter"),
       fill=c("red", "blue"))

i=1
while(i<=I){
set.seed(123)
this_pair=PAIR[i,]
this_l=which(BIN_FLAG==this_pair[1])
this_r=which(BIN_FLAG==this_pair[2])
this_l_vec=VEC[this_l,]
this_r_vec=VEC[this_r,]

start_point= this_l_vec[round(runif(1)*nrow(this_l_vec)),]
end_point= this_r_vec[round(runif(1)*nrow(this_r_vec)),]
points(start_point[1],start_point[2],col='red',pch=16,cex=2)
points(end_point[1],end_point[2],col='blue',pch=16,cex=2)
points(this_l_vec,col='black',pch=16,cex=0.3)
points(this_r_vec,col='black',pch=16,cex=0.3)

segments(start_point[1], start_point[2], end_point[1],end_point[2],col='grey40',lty=3,lwd=1)

#if(i!=1){
#text(x=max(VEC[,1]),y=max(VEC[,2]),label=paste0('TOP:',as.character(i-1)),pos=2,col='white')  
#text(x=max(VEC[,1]),y=min(VEC[,2]),label=paste0('SCORE:',as.character(SCORE[i-1])),pos=2,col='white')  }
#text(x=min(VEC[,1]),y=min(VEC[,2]),label=paste0('PERCENT:',as.character(round((i-1)/length(SNCMAT)*100)),'%'),pos=4,col='white') 
#}

i=i+1}
i=i-1
text(x=max(VEC[,1]),y=max(VEC[,2]),label=paste0('TOP:',as.character(i)),pos=2)  
text(x=max(VEC[,1]),y=min(VEC[,2]),label=paste0('SCORE:',as.character(SCORE[i])),pos=2) 
text(x=min(VEC[,1]),y=min(VEC[,2]),label=paste0('PERCENT:',as.character(round((i)/length(SNCMAT)*100)),'%'),pos=4)  
print(I)
I=I+1
dev.off()
}


 

















ALLOUT=c()


this_l_exp=PMAT[,83]
this_r_exp=PMAT[,30]

tag_list=c()
out_list=c()


i=1
while(i<=nrow(LR)){

this_l=LR[i,1]
this_r=LR[i,2]
this_tag=paste0(this_l,"_",this_r)
if(this_l %in% GENE & this_r %in% GENE){
    this_out=this_l_exp[which(names(this_l_exp)==this_l)]+this_r_exp[which(names(this_r_exp)==this_r)]
    tag_list=c(tag_list,this_tag)
    out_list=c(out_list,this_out)
    }  
  
i=i+1
}

names(out_list)=tag_list




sort_out_list=sort(out_list,decreasing=T)
sort_out_list[1:100]





library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))


BIN_FLAG=rep(NA,ncol(EXP))
i=1
while(i<=ncol(BIN)){
BIN_FLAG[BIN[,i]]=i
i=i+1
}


library('gplots')
heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))





##################
load('Seurat_EXP_cluster.Robj')
pbmc=EXP_cluster
pbmc@meta.data$bin=BIN_FLAG

TSNEPlot(object = pbmc,do.label=T,group.by='bin')

RS=apply(CMAT,2,sum)
LS=apply(CMAT,1,sum)


which(scale(RS)>1)
which(scale(LS)>1)

plot(RS,type='l')
points(LS,type='l',col='red')


which(RS==max(RS))
which(LS==max(LS))


LL=which(BIN_FLAG %in% c(59:83))
RR=which(BIN_FLAG %in% c(51,52,53,54,55,56,57,57,90,91,92,93))
pbmc@meta.data$lr=rep(NA,length(pbmc@ident))
pbmc@meta.data$lr[LL]='LL'
pbmc@meta.data$lr[RR]='RR'
pbmc@meta.data$ls=LS
TSNEPlot(object = pbmc,do.label=T,group.by='lr')










##########












































vis_gene='PDGFRA'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='SOX2'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='OLIG2'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='GFAP'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

vis_gene='CD83'
boxplot(as.numeric(EXP[which(GENE==vis_gene),])~BIN_FLAG)

plot(as.numeric(EXP[which(GENE==vis_gene),])~jitter(BIN_FLAG))

pbmc_data=EXP
pbmc= CreateSeuratObject(raw.data = pbmc_data, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc,do.plot=FALSE, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
PCNUM=50
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5)

PCUSE=1:50
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = pbmc,do.label=T)


pbmc@meta.data$bin=BIN_FLAG

TSNEPlot(object = pbmc,do.label=T,group.by='bin')


RS=apply(CMAT,2,sum)
LS=apply(CMAT,1,sum)

plot(RS,type='l')
points(LS,type='l',col='red')

which(scale(RS)>1)
which(scale(LS)>1)


which(RS==max(RS))
which(LS==max(LS))



which(CMAT[,20])

LL=which(BIN_FLAG %in% c(3))
LLL=which(BIN_FLAG %in% c(18))

RR=which(BIN_FLAG %in% c(20))
pbmc@meta.data$lr=rep(NA,length(pbmc@ident))
pbmc@meta.data$lr[LL]='LL'
pbmc@meta.data$lr[LLL]='LLL'
pbmc@meta.data$lr[RR]='RR'
pbmc@meta.data$ls=LS
TSNEPlot(object = pbmc,do.label=T,group.by='lr')


plot(BIN_FLAG,pbmc@dr$tsne@cell.embeddings[,1])
plot(BIN_FLAG,pbmc@dr$tsne@cell.embeddings[,2])

FeaturePlot(object = pbmc, cols.use = c("blue",'grey90', "red"), features.plot =c("PDGFRA",'BMI1','TP53','OLIG2','EGFR','GFAP'))

FeaturePlot(object = pbmc, cols.use = c("blue",'grey90', "red"), features.plot =c('ADAR','ADARB1','ADARB2'))
FeaturePlot(object = pbmc, cols.use = c("blue",'grey90', "red"), features.plot =c('ADAR','ADARB1','ADARB2'))
