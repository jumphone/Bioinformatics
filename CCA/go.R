
library(Seurat)
source('scRef.R')

tag_data=read.table('GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt',sep='\t',row.names=1,header=T)


pbmc.data=read.table('GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt',sep='\t',row.names=1,header=T)
#saveRDS(pbmc.data,'ALL.RDS')

MS=pbmc.data[,which(tag_data[,3]=='MS')]
CT=pbmc.data[,which(tag_data[,3]=='Ctrl')]

saveRDS(MS,'MS.RDS')
saveRDS(CT,'CT.RDS')

CPU=3
PCNUM=50
PCUSE=1:PCNUM

##########
MSD = CreateSeuratObject(raw.data = MS, min.cells = 0, min.genes = 0, project = "MS") 
MSD <- NormalizeData(object = MSD, normalization.method = "LogNormalize", scale.factor = 10000)
MSD <- ScaleData(object = MSD, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
MSD <- RunPCA(object = MSD, pcs.compute=PCNUM, pc.genes = rownames(MSD@data), do.print = FALSE)
MSD <- RunTSNE(object = MSD, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
MSX=MSD@dr$tsne@cell.embeddings
saveRDS(MSX,file='MSX.RDS')
##########
CTD = CreateSeuratObject(raw.data = CT, min.cells = 0, min.genes = 0, project = "CT") 
CTD <- NormalizeData(object = CTD, normalization.method = "LogNormalize", scale.factor = 10000)
CTD <- ScaleData(object = CTD, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
CTD <- RunPCA(object = CTD, pcs.compute=PCNUM, pc.genes = rownames(CTD@data), do.print = FALSE)
CTD <- RunTSNE(object = CTD, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
CTX=CTD@dr$tsne@cell.embeddings
saveRDS(CTX,file='CTX.RDS')


##########
##########
##########
##########
##########
##########
library(Seurat)
source('scRef.R')
library(pcaPP)

MS=readRDS("MS.RDS")
MSX=readRDS("MSX.RDS")
CT=readRDS("CT.RDS")
CTX=readRDS("CTX.RDS")


#GNUM=100
CNUM=100

.getGroup=function(X,TAG){
    DR=X
    RANK=rank(DR,ties.method='random')
    CUTOFF=CNUM #round(max(RANK)/GNUM)
    GROUP=rep('NA',length(RANK))
    i=1
    j=1
    while(i<=length(RANK)){
        GROUP[which(RANK==i)]=paste0(TAG,'_',as.character(j))
        if(i%%CUTOFF==1){j=j+1;print(j)}
        i=i+1}
    return(GROUP)
}



#####
MSG=.getGroup(MSX,'MS')
CTG=.getGroup(CTX,'CT')
########
MSR=.generate_ref(MS, cbind(MSG,MSG), min_cell=1) 
CTR=.generate_ref(CT, cbind(CTG,CTG), min_cell=1) 

COM=.simple_combine(MSR,CTR)
MSRC=COM$exp_sc_mat1
CTRC=COM$exp_sc_mat2

#saveRDS(MSR,file='MSR.RDS')
#saveRDS(CTR,file='CTR.RDS')
########

out = .get_cor( MSR, CTR, method='kendall',CPU=4, print_step=10)
tag1=.get_tag_max(out)
tag2=.get_tag_max(t(out))

V=c()
i=1
while(i<=nrow(tag1)){
    t1=tag1[i,1]
    t2=tag1[i,2]
    if(tag2[which(tag2[,1]==t2),2]==t1){V=c(V,i)}           
    i=i+1}
VP=tag1[V,]

C=c()
t=1
while(t<=nrow(VP)){
    this_c=out[which(rownames(out)==VP[t,2]),which(colnames(out)==VP[t,1])]
    C=c(C,this_c)
    t=t+1}
plot(C)
##########
VVP=VP[which(C>0.7),]
##########


TAG=c(MSG,CTG)
CON=c(rep('MS',length(MSG)),rep('CT',length(CTG)))
ALL=.simple_combine(MS,CT)$combine
pbmc=CreateSeuratObject(raw.data = ALL, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc@meta.data$tag=TAG
pbmc@meta.data$con=CON

MAP=rep('NA',length(TAG))
MAP[which(TAG %in% VVP[,1])]='MS'
MAP[which(TAG %in% VVP[,2])]='CT'
pbmc@meta.data$map=MAP


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"))

PCNUM=100
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)


VlnPlot(object = pbmc, features.plot = "PC1", group.by = "tag", do. return = F)
VlnPlot(object = pbmc, features.plot = "PC1", group.by = "con", do.return = F)




###############################

pdf('unMAP.pdf')
DimPlot(object =pbmc, reduction.use = "pca", group.by = "con",  pt.size = 0.5, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "pca", group.by = "map",  pt.size = 0.5, do.return = TRUE)
VlnPlot(object = pbmc, features.plot = "PC1", group.by = "map", do.return = TRUE)
VlnPlot(object = pbmc, features.plot = "PC2", group.by = "map", do.return = TRUE)
VlnPlot(object = pbmc, features.plot = "PC3", group.by = "map", do.return = TRUE)
dev.off()

#########

pbmc@dr$alnpca=pbmc@dr$pca
pbmc@dr$alnpca@key='APC'
#library(dtw)
#library(nloptr)
TAG=TAG
CON=CON
VALID_PAIR=VVP

index1=which(CON=='MS')
index2=which(CON=='CT')

ALL_COEF=c()
THIS_DR=1

while(THIS_DR<=ncol(pbmc@dr$pca@cell.embeddings)){
THIS_PC = pbmc@dr$pca@cell.embeddings[,THIS_DR]
M1=c()
M2=c()
S1=c()
S2=c()
i=1
while(i<=nrow(VALID_PAIR)){
    this_pair=VALID_PAIR[i,]
    this_index1=which(TAG %in% this_pair[1])
    this_index2=which(TAG %in% this_pair[2])
    this_m1=mean(THIS_PC[this_index1])
    this_m2=mean(THIS_PC[this_index2])
    this_s1=sd(THIS_PC[this_index1])
    this_s2=sd(THIS_PC[this_index2])
    M1=c(M1,this_m1)
    M2=c(M2,this_m2)
    i=i+1}
fit=lm(M2~M1)
this_coef=fit$coefficients
ALL_COEF=cbind(ALL_COEF,this_coef)
colnames(ALL_COEF)[THIS_DR]=as.character(THIS_DR)
#
pbmc@dr$alnpca@cell.embeddings[index1,THIS_DR]=ALL_COEF[1,THIS_DR]+pbmc@dr$pca@cell.embeddings[index1,THIS_DR]*ALL_COEF[2,THIS_DR]
#
THIS_DR=THIS_DR+1
}


pdf('MAP.pdf')
DimPlot(object =pbmc, reduction.use = "alnpca", group.by = "con",  pt.size = 0.5, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "alnpca", group.by = "map",  pt.size = 0.5, do.return = TRUE)
VlnPlot(object = pbmc, features.plot = "PC1", group.by = "map", do.return = TRUE)
VlnPlot(object = pbmc, features.plot = "PC2", group.by = "map", do.return = TRUE)
VlnPlot(object = pbmc, features.plot = "PC3", group.by = "map", do.return = TRUE)
VlnPlot(object = pbmc, features.plot = "APC50", group.by = "map", do.return = TRUE)
dev.off()




# VlnPlot(object = pbmc, features.plot = "PC1", group.by = "con", do.return = TRUE)
# PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2,group.by='con')

RATIO=0.8

PCUSE=which(ALL_COEF[2,] <1/RATIO & ALL_COEF[2,]>RATIO  )
pbmc <- RunTSNE(object = pbmc, reduction.use='alnpca',dims.use = PCUSE, do.fast = TRUE)

DimPlot(object =pbmc, reduction.use = "tsne", group.by = "con",  pt.size = 0.5, do.return = TRUE)

pdf('MAP.pdf')
DimPlot(object =pbmc, reduction.use = "tsne", group.by = "map",  pt.size = 0.5, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "tsne", group.by = "con",  pt.size = 0.5, do.return = TRUE)
dev.off()





























































































##################














mapped_index1 = which(TAG%in% VALID_PAIR[i,1])
mapped_index2 = which(TAG%in% VALID_PAIR[i,2])









X=this_coef

get_cost <- function(X){
     x0=X[1]
     x1=X[2]
     CHANGED_PC= x0+THIS_PC*x1
     COST=0
     i=1
     while(i<=nrow(VALID_PAIR)){
         this_pair=VALID_PAIR[i,]
         this_index1=which(TAG %in% this_pair[1])
         this_index2=which(TAG %in% this_pair[2])
         this_dist=(mean(CHANGED_PC[this_index1])-mean(THIS_PC[this_index2]))**2
         COST=COST+this_dist
     i=i+1
     }
     return(COST)
 }

get_cost(X)





optim(par = c(0.00001,-5),fn=get_cost,control = list(maxit=50000000))
 







this_pc_lst1 = THIS_PC[mapped_index1]
this_pc_lst2 = THIS_PC[mapped_index2]
this_pc_lst1 = sort(this_pc_lst1)
this_pc_lst2 = sort(this_pc_lst2)
this_aln = dtw(this_pc_lst1,this_pc_lst2,keep=TRUE)






























































######################################################################
######################################################################
######################################################################


library(dtw)

TAG=TAG
CON=CON
VALID_PAIR=VVP

index1=which(CON=='MS')
index2=which(CON=='CT')
mapped_index1 = which(TAG%in% VALID_PAIR[,1])
mapped_index2 = which(TAG%in% VALID_PAIR[,2])



THIS_DR=1
while(THIS_DR<=ncol(pbmc@dr$pca@cell.embeddings)){ 


THIS_PC = pbmc@dr$pca@cell.embeddings[,THIS_DR]
this_pc_lst1 = THIS_PC[mapped_index1]
this_pc_lst2 = THIS_PC[mapped_index2]
this_pc_lst1 = sort(this_pc_lst1)
this_pc_lst2 = sort(this_pc_lst2)
this_aln = dtw(this_pc_lst1,this_pc_lst2,keep=TRUE)


aln_pc_lst1 = rep(0, length(this_pc_lst1))
aln_pc_lst2 = rep(0, length(this_pc_lst2))

ii=1
while(ii<=length(aln_pc_lst1)){
    aln_pc_lst1[ii]=mean(which(this_aln$index1==ii))
    ii=ii+1}

ii=1
while(ii<=length(aln_pc_lst2)){
    aln_pc_lst2[ii]=mean(which(this_aln$index2==ii))
    ii=ii+1}


#plot(this_aln, type="threeway")
#plot(this_aln, type="twoway",offset=-2)

this_percent_lst1 <- ecdf(this_pc_lst1)
this_percent_lst2 <- ecdf(this_pc_lst2)

lst1_to_aln = function(X){Y=quantile(aln_pc_lst1, this_percent_lst1(X));return(Y)}
lst2_to_aln = function(X){Y=quantile(aln_pc_lst2, this_percent_lst2(X));return(Y)}

MAPPED_PC=THIS_PC
MAPPED_PC[index1]=lst1_to_aln(THIS_PC[index1])
MAPPED_PC[index2]=lst2_to_aln(THIS_PC[index2])

pbmc@dr$pca@cell.embeddings[,THIS_DR]=scale(MAPPED_PC)

print(THIS_DR)
THIS_DR=THIS_DR+1}





# VlnPlot(object = pbmc, features.plot = "PC1", group.by = "con", do.return = TRUE)
# PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2,group.by='con')
PCUSE=1:50
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)

DimPlot(object =pbmc, reduction.use = "tsne", group.by = "con",  pt.size = 0.5, do.return = TRUE)


MAP=rep('NA',length(TAG))
MAP[which(TAG %in% VVP[,1])]='MS'
MAP[which(TAG %in% VVP[,2])]='CT'
pbmc@meta.data$map=MAP

pdf('MAP.pdf')
DimPlot(object =pbmc, reduction.use = "tsne", group.by = "map",  pt.size = 0.5, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "tsne", group.by = "con",  pt.size = 0.5, do.return = TRUE)
dev.off()






















##########################################
##########################################
##########################################
##########################################
##########################################



X=c()
Y=c()
t=1
while(t<=nrow(VVP)){
    this_X=MSRC[,which(colnames(MSRC)==VVP[t,1])]
    this_Y=CTRC[,which(colnames(CTRC)==VVP[t,2])]
    X=cbind(X,this_X)
    Y=cbind(Y,this_Y)
    t=t+1}
colnames(X)=VVP[,1]
colnames(Y)=VVP[,2]



ALL=cbind(MS,CT)

PCA=princomp(ALL)


plot(PCA$loading[c(1:11),1], PCA$loading[c(12:22),1])








#install.packages('dtw')
library(dtw)

seq1=X[2,]
seq2=Y[2,]

seq1=seq1[order(seq1)]
seq2=seq2[order(seq2)]
A<-dtw(seq1,seq2,keep=TRUE);
plot(A, type="threeway")



plot(A, type="twoway",offset=-2);


























## A noisy sine wave as query
idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;

## A cosine is for template; sin and cos are offset by 25 samples
template<-cos(idx)

## Find the best match with the canonical recursion formula
library(dtw);
alignment<-dtw(query,template,keep=TRUE);
plot(alignment,type="threeway")


AL<-dtw(X[,1],Y[,1],keep=TRUE);
plot(AL,type="threeway")

plot(
    dtw(query,template,keep=TRUE,
        step=rabinerJuangStepPattern(6,"c")),
    type="twoway",offset=-2);












###############################
#install.packages('vegan')
library('vegan')

#install.packages('CCA')
library(CCA)


CC=.simple_combine(X,Y)$combine
PCA=princomp(CC)

COL=rep('black',22)
COL[c(1,12)]='red'
plot(PCA$loadings[,1],PCA$loading[,2],col=COL)

plot(PCA$loadings[1,],PCA$loadings[12,])



NX=X[,1:5]
NY=Y[,1:8]

XM=apply(NX,1,sum)
YM=apply(NY,1,sum)


NX=NX[which(XM>0 & YM>0),]
NY=NY[which(XM>0 & YM>0),]


IX=NX
IY=NY



#install.packages('candisc')
library(candisc)
ca=cancor(IX, IY)

U<-as.matrix(IX) %*% ca$coef$X
V<-as.matrix(IY) %*% ca$coef$Y





TMP=cca(IX,IY)
#cc1 <- cc(IX, IY)
#cc2 <- comput(IX, IY, cc1)
TMP$CCA$v

CC=.simple_combine(TMP$CCA$wa, IX)$combine



predict(TMP,newdata=IX)


#############
IY=t(NY)
IX=t(NX)
CCAX=cca(IY~IX)
CCAY=cca(IX~IY)
plot(CCAX)
plot(CCAY)



#plot(CCAX$CCA$u[,1], CCAY$CCA$u[,1])

#CCAX$CCA$u[1:3,1:3]

#test=data.frame(IX=t(NX[,c(1,2)]))
#tmp=predict(CCAX,newdata=test)


CC=.simple_combine(CCAX$CCA$v, NY)$combine










.getRankRatio=function(X){
    R=(rank(X,ties='min'))/max(rank(X,ties='min'))
    R[which(X==0)]=0
    return(R)
    }

################################################################

set.seed(123)
MIN_EXP_CUT=3
GROUP_SAMPLE_SIZE=100


OUT=matrix(ncol=2, nrow=GROUP_SAMPLE_SIZE* nrow(VVP))
i=1
while(i<=nrow(VVP)){
    in_this_group=VVP[i,1]
    out_this_group=VVP[i,2]
    in_this_group_exp=MS[,which(MSG %in% in_this_group)]
    out_this_group_exp=CT[,which(CTG %in% out_this_group)]
    j=1
    while(j<=GROUP_SAMPLE_SIZE){
        in_RR=apply(in_this_group_exp,1,sample,1)
        out_RR=apply(out_this_group_exp,1,sample,1)

        order_in_RR=in_RR[order(in_RR,decreasing=T)]
        order_out_RR=out_RR[order(out_RR,decreasing=T)]

        in_seq=names(order_in_RR)[which(order_in_RR>=MIN_EXP_CUT)]
        out_seq=names(order_out_RR)[which(order_out_RR>=MIN_EXP_CUT)]

        final_in_seq=paste(in_seq,collapse=',')
        final_out_seq=paste(out_seq,collapse=',')
        OUT[(i-1)*GROUP_SAMPLE_SIZE+j, 1]=final_in_seq
        OUT[(i-1)*GROUP_SAMPLE_SIZE+j, 2]=final_out_seq
        print(j)
        j=j+1}
    print(i)
    i=i+1}

colnames(OUT)=c('IN_SEQ','OUT_SEQ')
#write.table(OUT,file='TRAIN.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(OUT,file='TRAIN_noheader.txt',row.names=F,col.names=F,sep='\t',quote=F)
#write.table(OUT[,1],file='IN.txt',row.names=F,col.names=F,sep='\t',quote=F)
#write.table(OUT[,2],file='OUT.txt',row.names=F,col.names=F,sep='\t',quote=F)











#in_FF=in_this_group_exp[,1]
#out_FF=out_this_group_exp[,2]
#library(pcaPP)

#cor.fk(in_RR, out_RR)
#cor.fk(in_FF, out_FF)















X=c()
Y=c()
t=1
while(t<=nrow(VVP)){
    this_X=MSRC[,which(colnames(MSRC)==VVP[t,1])]
    this_Y=CTRC[,which(colnames(CTRC)==VVP[t,2])]
    X=cbind(X,this_X)
    Y=cbind(Y,this_Y)
    t=t+1}
colnames(X)=VVP[,1]
colnames(Y)=VVP[,2]


#install.packages('qualV')
library(qualV)

.getSeq=function(VALUE, NAME){
     O=order(-VALUE)
     N=NAME[O]
     VALUE=VALUE[O]
     OUT=c()
     i=1
     while(i<=length(VALUE) & VALUE[i]>0){
         OUT=c(OUT,N[i])
         i=i+1}
    return(OUT)
}

CC=X    
CC=CC*0

LL=c()
i=1
while(i<=ncol(X)){
    SX=.getSeq(X[,i],rownames(X))
    SY=.getSeq(Y[,i],rownames(Y))
    this_LL=LCS(SX,SY)$LCS
    CC[which(rownames(CC) %in% this_LL),i]=1
    print(i)
    i=i+1}

GG=apply(CC,1,sum)
used_gene=names(GG[which(GG>1)])


.getRankRatio=function(X){
    R=(rank(X,ties='min'))/max(rank(X,ties='min'))
    R[which(X==0)]=0
    return(R)
    }
RX=apply(X,2,.getRankRatio)
RY=apply(Y,2,.getRankRatio)

library(pcaPP)

plot(X[which(rownames(X) %in% used_gene),1], Y[which(rownames(X) %in% used_gene),1])
plot(RX[which(rownames(X) %in% used_gene),1], RY[which(rownames(X) %in% used_gene),1])

cor(RX[which(rownames(X) %in% used_gene),1], RY[which(rownames(X) %in% used_gene),1])
cor(RX[,1],RY[,1])

ULL=unique(LL)
LL=LCS(SX1,SY1)
LL$LCS[1]








plot(X[,1],Y[,1],pch=16)



################
.getRankRatio=function(X){
    R=(rank(X,ties='min'))/max(rank(X,ties='min'))
    R[which(X==0)]=0
    return(R)
    }
RX=apply(X,2,.getRankRatio)
RY=apply(Y,2,.getRankRatio)
RD=RY-RX




D=(Y-X)/CNUM

out=.get_log_p_sc_given_ref(MS, X, CPU=4, print_step=10)
tag=.get_tag_max(out)

table(tag[,2])

newMS=MS



































############
############
############


COM=.simple_combine(MS,CT)

matchedMS=COM$exp_sc_mat1
matchedCT=COM$exp_sc_mat2

mappedMS=matchedMS[,which(MSG %in% VVP[,1])]
mappedCT=matchedCT[,which(CTG %in% VVP[,2])]


rank_mappedMS=apply(mappedMS,2, .getRankRatio)
rank_mappedCT=apply(mappedMS,2, .getRankRatio)


rank_matchedMS=apply(matchedMS,2, .getRankRatio)
rank_matchedCT=apply(matchedCT,2, .getRankRatio)


#ecdf(c(0,rank_mappedMS[1,],1))(rank_matchedMS[1,1])
#quantile(rank_mappedCT,0.9)

#rank(rank_matchedMS[1,])


OUT=rank_matchedMS
OUT=OUT*0
i=1
while(i<=nrow(rank_matchedMS)){
    from_dist=rank_mappedMS[i,]
    to_dist=rank_mappedCT[i,]
    to_length=length(to_dist)
    to_sort=sort(to_dist)
    
    ECDF=ecdf(from_dist)
    j=1
    while(j<=ncol(rank_matchedMS)){
        this_p=ECDF(rank_matchedMS[i,j])
        this_out= to_sort[this_p*to_length]
        OUT[i,j]=this_out
        j=j+1}
    print(i)
    i=i+1}



.change_distribution <- function(exp_mat1, mapped_exp_mat1, mapped_exp_mat2, CPU=4, print_step=10){
    ##########
    library(parallel)
    #########
    print_step=print_step
    CPU=CPU
    exp_mat1=exp_mat1
    mapped_exp_mat1=mapped_exp_mat1
    mapped_exp_mat2=mapped_exp_mat2
    #######
    SINGLE <- function(i){
        from_dist=mapped_exp_mat1[i,]
        from_length=length(from_dist)
        to_dist=mapped_exp_mat2[i,]
        to_length=length(to_dist)
        to_sort=sort(to_dist)
        ECDF=ecdf(from_dist) 
        out=rep(0,from_length)
        j=1
        while(j<=ncol(exp_mat1)){
            this_exp=exp_mat1[i,j]
            if(this_exp!=0){
                this_p=ECDF(this_exp)
                this_out= to_sort[this_p*to_length]
            }else{this_out=0}
            out[j]=this_out
            j=j+1}        
        if(print_step==1){print(i)}else if(i%%print_step==1){print(i)}
        return(out)
        }
    ########
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:nrow(exp_mat1), SINGLE)
    stopCluster(cl)
    out = exp_mat1
    out = out*0
    i=1
    for(this_out in RUN){
        out[i,]=this_out
        i=i+1}   
    return(out)
}


exp_mat1=rank_matchedMS
mapped_exp_mat1=rank_mappedMS
mapped_exp_mat2=rank_mappedCT

out=.change_distribution(exp_mat1, mapped_exp_mat1, mapped_exp_mat2, CPU=4, print_step=10)




A=apply(out,2,.getRankRatio)
B=rank_matchedCT

MSA=CreateSeuratObject(raw.data = A, min.cells = 0, min.genes = 0, project = "MS")
CTB=CreateSeuratObject(raw.data = B, min.cells = 0, min.genes = 0, project = "MS")

pbmc=MergeSeurat(MSA,CTB,add.cell.id1='MS',add.cell.id2='CT')
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use =pbmc@var.genes)

PCNUM=50
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes)
PSUSE=1:50
pbmc <- RunTSNE(object = pbmc, dims.use =PSUSE, do.fast = TRUE)

TSNEPlot(object = pbmc,pt.size=0.1)
saveRDS(pbmc,'zfcombined.RDS')



write.table(pbmc@meta.data,file='zfcombined_meta.txt',row.names=T,col.names=T,sep='\t',quote=F)




pdf('zfcombined.pdf',width=10,height=7)
TSNEPlot(object = pbmc,pt.size=0.5)
dev.off()






