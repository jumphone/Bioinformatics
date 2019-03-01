
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

saveRDS(MSR,file='MSR.RDS')
saveRDS(CTR,file='CTR.RDS')
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

.getRankRatio=function(X){
    R=(rank(X,ties='min'))/max(rank(X,ties='min'))
    R[which(X==0)]=0
    return(R)
    }


COM=.simple_combine(MS,CT)

matchedMS=COM$exp_sc_mat1
matchedCT=COM$exp_sc_mat2

mappedMS=matchedMS[,which(MSG %in% VVP[,1])]
mappedCT=matchedCT[,which(CTG %in% VVP[,2])]


rank_mappedMS=apply(mappedMS,2, .getRankRatio)
rank_mappedCT=apply(mappedMS,2, .getRankRatio)


rank_matchedMS=apply(matchedMS,2, .getRankRatio)
rank_matchedCT=apply(matchedCT,2, .getRankRatio)



.change_distribution <- function(exp_mat1, mapped_exp_mat1, mapped_exp_mat2){}













################################################################
X=c()
Y=c()
t=1
while(t<=nrow(VVP)){
    this_X=MSRC[,which(colnames(MSRC)==VVP[t,1])]
    this_Y=CTRC[,which(colnames(CTRC)==VVP[t,2])]
    X=cbind(X,this_X)
    Y=cbind(Y,this_Y)
    t=t+1}


plot(X[,1],Y[,1],pch=16)





################
.getRankRatio=function(X){
    R=(rank(X,ties='min')-1)/max(rank(X,ties='min')-1)
    return(R)
    }

RX=apply(X,2,.getRankRatio)
RY=apply(Y,2,.getRankRatio)
RD=RY-RX

#D=dist(RD)
#H=hclust(D)




plot(RX[,1],RY[,1],pch=16)


library(pcaPP)


tmp=apply(RD,1,median)#*RX[,1]
tmpRX1=tmp+RX[,1]

tmpRX1


RtmpRX1=.getRankRatio(tmpRX1)
#plot(RtmpRX1,RX[,1])
cor.fk(RX[,1],RY[,1])
cor.fk(RtmpRX1,RY[,1])

par(mfrow=c(1,2))
#plot(RX[,1],RtmpRX1)
plot(RtmpRX1, RY[,1])
plot(RX[,1], RY[,1])










################
.getRankRatio=function(X){
    R=(rank(X,ties='min')-1)/max(rank(X,ties='min')-1)
    return(R)
    }

RX=apply(X,2,.getRankRatio)
RY=apply(Y,2,.getRankRatio)

X=RX
Y=RY

plot(X[,1],Y[,1])



D=Y-X
TD=t(D)

rownames(TD)=paste0(VVP[,1],'_',VVP[,2])


S=apply(TD,2,sd)
M=apply(TD,2,mean)

plot(S,M,pch=16)
plot(S,abs(M),pch=16)


XX=(X+M)




##########

.getRankRatio=function(X){
    R=(rank(X,ties='min')-1)/max(rank(X,ties='min')-1)
    return(R)
    }

RX=apply(X,2,.getRankRatio)
RY=apply(Y,2,.getRankRatio)

##########

X=RX
Y=RY

INTER_PCUT=0.05
COEF_PCUT=0.05

INTER=c()
COEF=c()

i=1
while(i<=nrow(X)){
    fit=lm(Y[i,]~X[i,])
    this_sum=summary(fit)   
    this_inter=fit$coefficients[1]
    this_coef=fit$coefficients[2]
    this_inter_p=1
    this_coef_p=1
    if(!is.na(this_inter) & this_inter!=0){
        this_inter_p=this_sum$coefficients[1,4]}
    if(!is.na(this_coef) & this_coef!=0 ){
        this_coef_p=this_sum$coefficients[2,4]}
    out_inter=0
    out_coef=1
    if((!is.na(this_inter)) & (this_inter_p<INTER_PCUT) ){out_inter=this_inter}
    if((!is.na(this_coef)) & (this_coef_p<COEF_PCUT) ){out_coef=this_coef}
    INTER=c(INTER, out_inter)
    COEF=c(COEF, out_coef)
    print(i)
    i=i+1}




XX=X*COEF+INTER
library('pcaPP')
i=3
cor.fk(X[,i],Y[,i])
cor.fk(XX[,i],Y[,i])






###################
library('pcaPP')


.getRankRatio=function(X){
    R=(rank(X,ties='min')-1)/max(rank(X,ties='min')-1)
    return(R)
    }

RX=apply(X,2,.getRankRatio)
RY=apply(Y,2,.getRankRatio)




RR=apply(cbind(tmpX,tmpY),1,max)


tmpC=cor.fk(tmpX,tmpY)

deltaC=c()
i=1
while(i<=length(tmpX)){
    thisX=tmpX[i]
    thisY=tmpY[i]
    LX=tmpX
    LY=tmpY
    LX[i]=thisY
    LY[i]=thisX
    thisC1=cor.fk(LX,tmpY)
    thisC2=cor.fk(LY,tmpX)
    thisdC=(thisC1+thisC2)/2-tmpC
    deltaC=c(deltaC,thisdC)
    print(i)
    i=i+1}
 

















