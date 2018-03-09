ALLCNV=as.matrix(read.table('relative_ALLCNV.txt',row.names=1,header=T,sep='\t' ))
CELL_COLOR=as.character(read.table('CELL_COLOR.txt',header=F,sep='\t' )[,1])
CHROM_COLOR=as.character(read.table('CHROM_COLOR.txt',header=F,sep='\t' )[,1])
CHROM_TAG=as.factor(read.table('CHROM_TAG.txt',header=F,sep='\t' )[,1])



ALLCNV[which(ALLCNV >3 )]=3
ALLCNV[which(ALLCNV < -3 )]= -3

pdf('CNV.pdf',width=40,height=40)

par(mfrow=c(length(ALLCNV[,1]),1))

i=1
while(i<=length(ALLCNV[,1])){

boxplot(as.numeric(ALLCNV[i,])~CHROM_TAG,ylim=c(-3,3),pch=16,main=rownames(ALLCNV)[i])
	i=i+1
}

par(mfrow=c(length(ALLCNV[,1]),1))

i=1
while(i<=length(ALLCNV[,1])){

plot(as.numeric(ALLCNV[i,]),col=CHROM_COLOR,ylim=c(-3,3),pch=16,main=rownames(ALLCNV)[i])
	i=i+1
}
#library(gplots)
#heatmap.2( ALLCNV,scale = c("none"),col=colorRampPalette(c("blue", "white","red"))(100),
#       Rowv=F,Colv=F,dendrogram='none',trace='none',labCol=F,labRow=F, RowSideColors=CELL_COLOR, ColSideColors=CHROM_COLOR)

#heatmap.2( ALLCNV,scale = c("none"),col=colorRampPalette(c("blue", "white","red"))(100),
#       Rowv=T,Colv=F,dendrogram='none',trace='none',labCol=F,labRow=F, RowSideColors=CELL_COLOR, ColSideColors=CHROM_COLOR)


dev.off()

