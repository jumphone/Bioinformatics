library('gplots')

rawdata=read.delim('combined.txt',header=T,row.names=1,sep='\t')
EXP=rawdata[,c(1,2,3,4)]
VAR=apply(EXP,1,var)
#VAR_BIG10=which(log(VAR,10) >=0)
MEAN=apply(EXP,1,mean)
analyzed_gene=which(MEAN >=10 & VAR>=10)


FINAL_EXP=log(EXP[analyzed_gene,]+1,2)
scale_center=function(x){return(scale(x,center=T,scale=F))}
LIMIT=1
S_FINAL_EXP=t(apply(FINAL_EXP,1,scale_center))
colnames(S_FINAL_EXP)=colnames(FINAL_EXP)
S_FINAL_EXP[which(S_FINAL_EXP >=LIMIT)]=LIMIT
S_FINAL_EXP[which(S_FINAL_EXP <=-LIMIT)]=-LIMIT

#hmr=heatmap.2(as.matrix(S_FINAL_EXP),scale=c("none"),dendrogram='row',trace='none',col=colorRampPalette(c('blue','white','red')))
#hmr=heatmap.2(as.matrix(FINAL_EXP),scale=c("row"),dendrogram='row',trace='none',col=colorRampPalette(c('blue','white','red')))
DIST=dist(S_FINAL_EXP)
HCLUST=hclust(DIST)
N=cutree(HCLUST,k=7)

COL=rep('red',length(S_FINAL_EXP[,1]))
COL[which(N==1)]='red'
COL[which(N==2)]='blue'
COL[which(N==3)]='green'
COL[which(N==4)]='yellow'
COL[which(N==5)]='purple'
COL[which(N==6)]='black'
COL[which(N==7)]='pink'

pdf('HEATMAP.pdf',width=30,height=30)
hmr=heatmap.2(as.matrix(S_FINAL_EXP),scale=c("none"),dendrogram='row',trace='none',col=colorRampPalette(c('blue','white','red')),RowSideColors=COL)

boxplot(S_FINAL_EXP[which(N==1),],main='1',col=c('red','red','blue','blue','green','green','purple','purple'))
boxplot(S_FINAL_EXP[which(N==2),],main='2',col=c('red','red','blue','blue','green','green','purple','purple'))
boxplot(S_FINAL_EXP[which(N==3),],main='3',col=c('red','red','blue','blue','green','green','purple','purple'))
boxplot(S_FINAL_EXP[which(N==4),],main='4',col=c('red','red','blue','blue','green','green','purple','purple'))
boxplot(S_FINAL_EXP[which(N==5),],main='5',col=c('red','red','blue','blue','green','green','purple','purple'))
boxplot(S_FINAL_EXP[which(N==6),],main='6',col=c('red','red','blue','blue','green','green','purple','purple'))
boxplot(S_FINAL_EXP[which(N==7),],main='7',col=c('red','red','blue','blue','green','green','purple','purple'))



write.table(EXP[ rownames(EXP) %in% rownames(S_FINAL_EXP)[which(N==1)],     ],'1.txt',quote=F,sep='\t')
write.table(EXP[ rownames(EXP) %in% rownames(S_FINAL_EXP)[which(N==2)],     ],'2.txt',quote=F,sep='\t')
write.table(EXP[ rownames(EXP) %in% rownames(S_FINAL_EXP)[which(N==3)],     ],'3.txt',quote=F,sep='\t')
write.table(EXP[ rownames(EXP) %in% rownames(S_FINAL_EXP)[which(N==4)],     ],'4.txt',quote=F,sep='\t')
write.table(EXP[ rownames(EXP) %in% rownames(S_FINAL_EXP)[which(N==5)],     ],'5.txt',quote=F,sep='\t')
write.table(EXP[ rownames(EXP) %in% rownames(S_FINAL_EXP)[which(N==6)],     ],'6.txt',quote=F,sep='\t')
write.table(EXP[ rownames(EXP) %in% rownames(S_FINAL_EXP)[which(N==7)],     ],'7.txt',quote=F,sep='\t')



dev.off()


