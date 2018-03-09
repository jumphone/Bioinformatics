ALLCNV=as.matrix(read.table('relative_ALLCNV.txt',row.names=1,header=T,sep='\t' ))
CELL_COLOR=as.character(read.table('CELL_COLOR.txt',header=F,sep='\t' )[,1])
CHROM_COLOR=as.character(read.table('CHROM_COLOR.txt',header=F,sep='\t' )[,1])



ALLCNV[which(ALLCNV >3 )]=3
ALLCNV[which(ALLCNV < -3 )]= -3

pdf('HEATMAP.pdf',width=30,height=7)


library(gplots)
heatmap.2( ALLCNV,scale = c("none"),col=colorRampPalette(c("blue", "white","red"))(100),
       Rowv=F,Colv=F,dendrogram='none',trace='none',labCol=F,labRow=F, RowSideColors=CELL_COLOR, ColSideColors=CHROM_COLOR)

heatmap.2( ALLCNV,scale = c("none"),col=colorRampPalette(c("blue", "white","red"))(100),
       Rowv=T,Colv=F,dendrogram='none',trace='none',labCol=F,labRow=F, RowSideColors=CELL_COLOR, ColSideColors=CHROM_COLOR)


dev.off()

