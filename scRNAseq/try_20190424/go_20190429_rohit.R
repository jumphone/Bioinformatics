library(Seurat)
library(dplyr)
library(Matrix)
#pbmc.data <- read.table("MGH54_mat.txt", sep='\t',row.names=1, header=T)

#load('Seurat_EXP_cluster.Robj')
#pbmc.raw.data=as.matrix(EXP_cluster@raw.data[,which(colnames(EXP_cluster@raw.data) %in% colnames(EXP_cluster@scale.data))])
#pbmc.data=as.matrix(EXP_cluster@scale.data)


load('NewWeinanWorkspace.RData')

pbmc=mmap
pbmc.raw.data=as.matrix(pbmc@raw.data[,which(colnames(pbmc@raw.data) %in% colnames(pbmc@scale.data))])
pbmc.data=as.matrix(pbmc@scale.data)

used=c(1:ncol(pbmc.raw.data))

###########
#pbmc.nonscale.data=as.matrix(EXP_cluster@data)[,used]
#pbmc.data=tmp@scale.data #pbmc.data[which(rownames(pbmc.data) %in% tmp@var.genes),]

#pbmc.data=pbmc.data[which(rownames(pbmc.data) %in% tmp@var.genes),]
#LR=read.table('RL.txt',header=T,sep='\t')



source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#rm(EXP_cluster)
#gc()

ONE=.data2one(pbmc.raw.data, rownames(pbmc.data), CPU=4, PCNUM=50, SEED=123,  PP=30)

WINDOW= round(length(ONE)/200)
#WINDOW=300
RANK=rank(ONE)
LENGTH=length(ONE)
BIN=c()
i=1
while(i<=LENGTH/WINDOW  ){
this_index=which((i-1)*WINDOW< RANK & i*WINDOW>=RANK)
BIN=cbind(BIN,this_index)
i=i+1
}

pbmc=mmap
VEC=pbmc@dr$tsne@cell.embeddings

BIN_FLAG=rep(NA,length(pbmc@ident))
i=1
while(i<=ncol(BIN)){
BIN_FLAG[used][BIN[,i]]=i
i=i+1
}

pbmc@meta.data$bin=BIN_FLAG
DimPlot(pbmc,group.by='bin',reduction.use='tsne',do.label=T)

saveRDS(BIN,'BIN.RDS')
saveRDS(ONE,'ONE.RDS')
################################################

LR=read.table('RL_mouse.txt',header=T,sep='\t')

EXP=pbmc.data
GENE=rownames(EXP)
ALL=c(as.character(LR[,1]),as.character(LR[,2]) )

permu_gene_index=which(GENE %in% ALL)

set.seed(123)
TIME=10000
MEAN=  matrix(nrow=nrow(EXP[permu_gene_index,]),ncol=TIME)#EXP[permu_gene_index,c(1:TIME)]*0
MEAN[which(is.na(MEAN))]=0
rownames(MEAN)=rownames(EXP[permu_gene_index,])
colnames(MEAN)=as.character(c(1:TIME))
i=1
while(i<=TIME){
this_index=sample(c(1:ncol(EXP)),WINDOW)
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
#CMAT=readRDS('CMAT.RDS')
#BIN=readRDS('BIN.RDS')
#library(Seurat)
#load('Seurat_EXP_cluster.Robj')

NCMAT=as.numeric(CMAT)
SNCMAT=sort(NCMAT,decreasing=T)
length(SNCMAT)*0.05

#TOP=round(length(SNCMAT)*0.05)
TOP=1000
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
#pbmc@meta.data$bin=BIN_FLAG
VEC=pbmc@dr$tsne@cell.embeddings

BIN_FLAG=rep(NA,ncol(as.matrix(pbmc@data)))
i=1
while(i<=ncol(BIN)){
BIN_FLAG[used][BIN[,i]]=i
i=i+1
}

pbmc@meta.data$bin=BIN_FLAG
#library('gplots')
#heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
#DimPlot(pbmc,group.by='bin',reduction.use='tsne')


'''
par(mfrow=c(2,1))
hist(as.numeric(CMAT),xlab='SCORE',breaks=100,freq=F,main='Histogram of SCORE')
abline(v=CUTOFF,col='red')
text(x=CUTOFF,y=0,pos=4,col='red',label=as.character(round(CUTOFF)))
'''


png('COMMUNICATION_TUMOR.png',width=1200,height=1000)
plot(VEC,col='grey80',pch=16,cex=0.3,main=paste0('TOP:',as.character(TOP), 
                                                 '; PERCENT:', as.character(round(TOP/length(SNCMAT)*100)),'%',
                                                 '; CUTOFF:',as.character(round(CUTOFF))))

legend("topleft", legend=c("Ligand", "Recepter"),
       fill=c("yellowgreen", "cornflowerblue"))

SIZE=c()
i=1
while(i<=nrow(PAIR)){
#set.seed(123)
this_pair=PAIR[i,]
this_l=which(BIN_FLAG==this_pair[1])
this_r=which(BIN_FLAG==this_pair[2])
this_l_vec=VEC[this_l,]
this_r_vec=VEC[this_r,]

#used_index=0.5
#start_point= this_l_vec[round(used_index*nrow(this_l_vec)),]
#end_point= this_r_vec[round(used_index*nrow(this_r_vec)),]
library(cluster)
  
start_point=pam(this_l_vec, 1)$medoids
end_point= pam(this_r_vec, 1)$medoids
#transparent_ratio =  
#points(start_point[1],start_point[2],pch=16,cex=2,col=rgb(255, 0, 0, transparent_ratio, maxColorValue=255) )
#points(end_point[1],end_point[2],pch=16,cex=2,col=rgb(0, 0, 255, transparent_ratio, maxColorValue=255) )
size_ratio = (nrow(PAIR)-i+1)/nrow(PAIR)
SIZE=c(SIZE,size_ratio)
base_size=4
points(start_point[1],start_point[2],pch=16,cex=base_size*size_ratio, col='yellowgreen' )
points(end_point[1],end_point[2],pch=16,cex=base_size*size_ratio, col='cornflowerblue' )
 
points(this_l_vec,col='grey50',pch=16,cex=0.3)
points(this_r_vec,col='grey50',pch=16,cex=0.3)
#text(x= this_l_vec[,1],y=this_l_vec[,2],label= as.character(this_pair[1]),cex=0.3,col='grey50')
#text(x= this_r_vec[,1],y=this_r_vec[,2],label= as.character(this_pair[2]),cex=0.3,col='grey50')
  
#coral1
text_col='red'
text_cex=1
text(x=start_point[1],y=start_point[2],label=as.character(this_pair[1]),pos=as.numeric(this_pair[1])%%4+1, col=text_col,cex=text_cex)  
text(x=end_point[1],y=end_point[2],label=as.character(this_pair[2]),pos=as.numeric(this_pair[2])%%4+1, col=text_col,cex=text_cex)  
  

segments(start_point[1], start_point[2], end_point[1],end_point[2],col='grey40',lty=3,lwd=1)
i=i+1}
dev.off()









DimPlot(pbmc, reduction.use='tsne', group.by='ident', pt.size=0.1,do.label=T)

TAG=pbmc@meta.data$RohitAnnotation
TAG[which(TAG=='Border\nMacrophages')]='Border_Macrophages'
pbmc@meta.data$newtag=TAG
#DimPlot(pbmc, reduction.use='tsne', group.by='newtag', pt.size=0.1,do.label=T)


TT=table(TAG,BIN_FLAG)
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
tag=.get_tag_max(TT)


PP=c()
i=1
while(i<=nrow(PAIR)){
this_pair=PAIR[i,]
this_tag1=tag[which(tag[,1]==as.character(this_pair[1])),2]
this_tag2=tag[which(tag[,1]==as.character(this_pair[2])),2]
PP=cbind(PP,c(this_tag1,this_tag2))
i=i+1}
PP=t(PP)

OUTPUT=cbind(PAIR,PP,SIZE)
colnames(OUTPUT)=c('L','R','LT','RT','SIZE')
write.table(OUTPUT,file='OUTPUT.txt',row.names=F,col.names=T,sep='\t',quote=F)
saveRDS(OUTPUT,file='OUTPUT.RDS')

########################

png('LABEL.png',width=1200,height=1000)
DimPlot(pbmc,group.by='RohitAnnotation',reduction.use='tsne',do.label=T)
dev.off()



png('PIE.png',width=1200,height=1200)
TOT=paste0(OUTPUT[,3],'_to_',OUTPUT[,4])
par(mar=c(15,15,15,15))
pie(sort(table(TOT)))
dev.off()









VP=OUTPUT[which(OUTPUT[,3]=='Microglia' & OUTPUT[,4]=='Tumor Cells'),]

this_l_exp=apply(PMAT[,as.numeric(VP[,1])],1,mean)
this_r_exp=apply(PMAT[,as.numeric(VP[,2])],1,mean)

tag_list=c()
out_list=c()
l_list=c()
r_list=c()

i=1
while(i<=nrow(LR)){

this_l=LR[i,1]
this_r=LR[i,2]
this_tag=paste0(this_l,"_",this_r)
if(this_l %in% GENE & this_r %in% GENE){
    this_out=this_l_exp[which(names(this_l_exp)==this_l)]+this_r_exp[which(names(this_r_exp)==this_r)]
    this_l_out=this_l_exp[which(names(this_l_exp)==this_l)]
    this_r_out=this_r_exp[which(names(this_r_exp)==this_r)]
    tag_list=c(tag_list,this_tag)
    out_list=c(out_list,this_out)
    l_list=c(l_list,this_l_out)
    r_list=c(r_list,this_r_out)
    }  
  
i=i+1
}

names(out_list)=tag_list

png('Microglia_Tumor.png',width=1200,height=1000)
set.seed(1234567)
names(l_list)=paste0(tag_list,'_L')
names(r_list)=paste0(tag_list,'_R')

#plot(r_list,l_list,pch=16,xlab='EXP of Receptor in Tumor Cell',ylab='EXP of Ligand in Microglia',xlim=c(-0.3,4),ylim=c(-0.3,4))
plot(r_list,l_list,pch=16,xlab='EXP of Receptor in Tumor Cell',ylab='EXP of Ligand in Microglia')
text(r_list,l_list,label=tag_list,pos=sample(c(1,2,3,4),length(l_list),replace = TRUE))
dev.off()





VP=OUTPUT[which(OUTPUT[,3]=='Tumor Cells' & OUTPUT[,4]=='Monocyte-Derived Cells'),]

this_l_exp=apply(PMAT[,as.numeric(VP[,1])],1,mean)
this_r_exp=apply(PMAT[,as.numeric(VP[,2])],1,mean)

tag_list=c()
out_list=c()
l_list=c()
r_list=c()

i=1
while(i<=nrow(LR)){

this_l=LR[i,1]
this_r=LR[i,2]
this_tag=paste0(this_l,"_",this_r)
if(this_l %in% GENE & this_r %in% GENE){
    this_out=this_l_exp[which(names(this_l_exp)==this_l)]+this_r_exp[which(names(this_r_exp)==this_r)]
    this_l_out=this_l_exp[which(names(this_l_exp)==this_l)]
    this_r_out=this_r_exp[which(names(this_r_exp)==this_r)]
    tag_list=c(tag_list,this_tag)
    out_list=c(out_list,this_out)
    l_list=c(l_list,this_l_out)
    r_list=c(r_list,this_r_out)
    }  
  
i=i+1
}

names(out_list)=tag_list

png('Tumor_MonocyteDerived.png',width=1200,height=1000)
#set.seed(12345)
names(l_list)=paste0(tag_list,'_L')
names(r_list)=paste0(tag_list,'_R')

#plot(r_list,l_list,pch=16,xlab='EXP of Receptor in Tumor Cell',ylab='EXP of Ligand in Microglia',xlim=c(-0.3,4),ylim=c(-0.3,4))
plot(r_list,l_list,pch=16,xlab='EXP of Receptor in MonocyteDerived',ylab='EXP of Ligand in Tumor')
text(r_list,l_list,label=tag_list,pos=sample(c(1,2,3,4),length(l_list),replace = TRUE))
dev.off()





















sort_out_list=sort(out_list,decreasing=T)
sort_out_list[1:100]

write.table(file='SIG_PAIR.txt',sort_out_list,col.names=F,row.names=T,quote=F,sep='\t')



pdf('EXP.pdf',width=16,height=6)
vis_gene='Thbs1'
par(mfrow=c(1,2))
boxplot(as.numeric(EXP[which(rownames(EXP)==vis_gene),])~BIN_FLAG[used],ylab='Nomalized_Expression',main=vis_gene,pch=16,las=2)
boxplot(as.numeric(pbmc.raw.data[which(rownames(pbmc.raw.data)==vis_gene),])~BIN_FLAG[used],ylab='Raw_Expression',main=vis_gene,pch=16,las=2)

vis_gene='Sdc4'
par(mfrow=c(1,2))
boxplot(as.numeric(EXP[which(rownames(EXP)==vis_gene),])~BIN_FLAG[used],ylab='Normalized_Expression',main=vis_gene,pch=16,las=2)
boxplot(as.numeric(pbmc.raw.data[which(rownames(pbmc.raw.data)==vis_gene),])~BIN_FLAG[used],ylab='Raw_Expression',main=vis_gene,pch=16,las=2)

vis_gene='Serpine1'
par(mfrow=c(1,2))
boxplot(as.numeric(EXP[which(rownames(EXP)==vis_gene),])~BIN_FLAG[used],ylab='Normalized_Expression',main=vis_gene,pch=16,las=2)
boxplot(as.numeric(pbmc.raw.data[which(rownames(pbmc.raw.data)==vis_gene),])~BIN_FLAG[used],ylab='Raw_Expression',main=vis_gene,pch=16,las=2)

vis_gene='Lrp1'
par(mfrow=c(1,2))
boxplot(as.numeric(EXP[which(rownames(EXP)==vis_gene),])~BIN_FLAG[used],ylab='Normalized_Expression',main=vis_gene,pch=16,las=2)
boxplot(as.numeric(pbmc.raw.data[which(rownames(pbmc.raw.data)==vis_gene),])~BIN_FLAG[used],ylab='Raw_Expression',main=vis_gene,pch=16,las=2)

vis_gene='Lrp2'
par(mfrow=c(1,2))
boxplot(as.numeric(EXP[which(rownames(EXP)==vis_gene),])~BIN_FLAG[used],ylab='Normalized_Expression',main=vis_gene,pch=16,las=2)
boxplot(as.numeric(pbmc.raw.data[which(rownames(pbmc.raw.data)==vis_gene),])~BIN_FLAG[used],ylab='Raw_Expression',main=vis_gene,pch=16,las=2)

dev.off()












#plot(as.numeric(pbmc.raw.data[which(GENE==vis_gene),])~jitter(BIN_FLAG[used]))
