source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')


library(Seurat)
library(dplyr)
library(Matrix)
load('HDN_Wkspace.RData')

pbmc=tumor2
pbmc.raw.data=as.matrix(pbmc@raw.data[,which(colnames(pbmc@raw.data) %in% colnames(pbmc@scale.data))])
pbmc.data=as.matrix(pbmc@scale.data)

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
ONE=.data2one(pbmc.raw.data, rownames(pbmc.data), CPU=4, PCNUM=50, SEED=123,  PP=30)
#ONE=readRDS('ONE.RDS')


OUT=getBIN(ONE)
BIN=OUT$BIN
BINTAG=OUT$TAG

pbmc@meta.data$bin=BINTAG
png('ID.png',width=1200,height=1000)
DimPlot(pbmc,group.by='bin',reduction.use='umap',do.label=T)
dev.off()



EXP=pbmc.data
LR=read.table('RL_mouse.txt',header=T,sep='\t')

MEAN=getMEAN(EXP, LR)

PMAT=getPMAT(EXP, LR, BIN, MEAN)
CMAT=getCMAT(EXP,LR,PMAT)



OUT=getPAIR(CMAT)

PAIR=OUT$PAIR
SCORE=OUT$SCORE
RANK=OUT$RANK







TAG=as.character(pbmc@ident)#pbmc@meta.data$RohitAnnotation
TAG[which(TAG=='Pericyte/\nFibroblast')]='Pericyte Fibroblast'

TT=table(TAG,BINTAG)
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


YES=cbind(PAIR,PP)
TOT=paste0(YES[,3],'_to_',YES[,4])
CCLR=names(sort(table(TOT),decreasing=T)[which(sort(table(TOT),decreasing=T)>3)])


par(mfrow=c(round((length(CCLR)+1)/2),2))
i=1
while(i<=length(CCLR)){
this_cclr=CCLR[i]
this_p=ks.test(which(TOT==this_cclr),1:length(TOT),alternative='greater')$p.value
this_p=signif(this_p, digits = 2)
this_p=format(this_p, scientific = T)
plot(main=paste0(as.character(i),': ',this_cclr,'; KS p-value=',this_p),x=which(TOT==this_cclr),y=rep(1,length(which(TOT==this_cclr))),type='h',ylim=c(0,1),xlim=c(0,length(TOT)),xlab='RANK',ylab='',col=i)
 
i=i+1}



####################################
VEC=pbmc@dr$umap@cell.embeddings
PAIR=PAIR[1:1000,]
BIN_FLAG=BINTAG
png('COMMUNICATION.png',width=1200,height=1000)
plot(VEC,col='grey80',pch=16,cex=0.3,main=paste0('TOP:',as.character(TOP), 
                                                 '; PERCENT:', as.character(round(TOP/length(SNCMAT)*100)),'%',
                                                 '; CUTOFF:',as.character(round(CUTOFF))))

legend("topleft", legend=c("Ligand", "Recepter"),
       fill=c("green", "blue"))

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
size_ratio = (nrow(PAIR)-i+1)/nrow(PAIR)
SIZE=c(SIZE,size_ratio)
base_size=4
#points(start_point[1],start_point[2],pch=16,cex=base_size*size_ratio, col='yellowgreen' )
#points(end_point[1],end_point[2],pch=16,cex=base_size*size_ratio, col='cornflowerblue' )
transparent_ratio =150
points(start_point[1],start_point[2],pch=16,cex=base_size*size_ratio, col=rgb(0, 255, 0, transparent_ratio, maxColorValue=255)  )
points(end_point[1],end_point[2],pch=16,cex=base_size*size_ratio, col=rgb(0, 0, 255, transparent_ratio, maxColorValue=255) )

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





