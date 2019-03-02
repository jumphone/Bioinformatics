
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






