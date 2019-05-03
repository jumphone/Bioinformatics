source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
D4=read.table('DIPG4_1_bismark_bt2_pe.bedGraph.gene.txt',sep='\t',row.names=1,header=F)
D13=read.table('DIPG13_1_bismark_bt2_pe.bedGraph.gene.txt',sep='\t',row.names=1,header=F)
D1=read.table('DIPGC1_1_bismark_bt2_pe.bedGraph.gene.txt',sep='\t',row.names=1,header=F)
D2=read.table('DIPGC2_1_bismark_bt2_pe.bedGraph.gene.txt',sep='\t',row.names=1,header=F)


D1=as.matrix(D1)
D2=as.matrix(D2)
D4=as.matrix(D4)
D13=as.matrix(D13)

C1=.simple_combine(cbind(D1,D1),cbind(D2,D2))$combine[,c(1,3)]
C2=.simple_combine(cbind(D4,D4),cbind(D13,D13))$combine[,c(1,3)]
C3=.simple_combine(C1,C2)$combine
colnames(C3)=c('DIPGC1','DIPGC2','DIPG4','DIPG13')
write.table(C3,'combine.txt',sep='\t',quote=F,row.names=T,col.names=T)
colnames(C3)=c('DIPGC1','DIPGC2','DIPG4','DIPG13')
C3=as.matrix(C3)








library(Seurat)

EXP = CreateSeuratObject(raw.data =C3, min.cells = 0, min.genes=0)
allgene=rownames(EXP@data)
EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)
EXP@meta.data$name=colnames(C3)


EXP = ScaleData(object = EXP,  genes.use = allgene)
heatmap(cor(C3))
PCNUM=3
EXP <- RunPCA(object = EXP, pc.genes = allgene, do.print = TRUE, pcs.print = 1:3,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )

pdf('PCA.pdf',width=5,height=4)
PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2, do.label=F, no.legend =F, do.hover = F,group.by='name',pt.size=5)
PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2, do.label=T, no.legend =F, do.hover = F,group.by='name')
dev.off()
#OUT=PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2, do.label=F, no.legend =T, do.hover = T)



OUT=as.matrix(EXP@data)
#VAR=apply(OUT,1,var)
#used_index=which(VAR>1)#which(VAR>=sort(VAR,decreasing=T)[2000])


pdf('HEAT.pdf',width=5,height=6)
library('gplots')
#heatmap.2(OUT[used_index,],scale=c("row"),dendrogram='column',Colv=T,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,5),labRow='')
heatmap.2(cor(OUT,method='spearman'),scale=c("none"),dendrogram='column',Colv=T,trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,5),labRow='')
dev.off()
