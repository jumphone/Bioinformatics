SUR=read.table('SUR.txt',sep='\t',header=T)
used=which(!is.na(SUR[,2]))

library(survival)
library(survminer)

score=SUR[used,2]
this_sur=SUR[used,3]

plot(this_sur,score, xlab='OS (month)', ylab='Reponse Rate (%)', pch=16)
abline(h=0,col='red')
cor.test(this_sur, score, method='spearman')

highcut=0
lowcut=0

TYPE=rep('MED',length(used))
TYPE[which(score> highcut)]='High'
TYPE[which(score< lowcut)]='Low'
surtime=this_sur[which(TYPE!='MED')]
surevent=rep(1,length(which(TYPE!='MED')))
surtype=TYPE[which(TYPE!='MED')]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
ggsurvplot(fit, pval = TRUE)
surv_pvalue(fit)



####################################

mut_data=read.table('COM.txt',header=T,sep='\t',row.names=1)
library(Seurat)
pbmc <- CreateSeuratObject(raw.data = mut_data, min.cells = 0, min.genes = 0, project = "PBTR")
pbmc@meta.data$nMut=pbmc@meta.data$nGene
VlnPlot(object = pbmc, features.plot = c("nMut"), nCol = 1)
dim(pbmc@data)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc,do.plot=F, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
#50,156
pbmc <- ScaleData(object = pbmc)
PCNUM=20
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:7, do.fast = TRUE,perplexity=5)
TSNEPlot(object = pbmc,pt.size=3,group.by='orig.ident')

TSNE_VEC=pbmc@dr$tsne@cell.embeddings
D=dist(TSNE_VEC)
H=hclust(D)
C=cutree(H, k=4) 
pbmc@meta.data$C=C
TSNEPlot(object = pbmc,pt.size=3 , do.label=T, group.by ='C')
save.image('PBTR_Seurat.Robj')

########

HEATDATA=t(table(pbmc@meta.data$C,names(pbmc@ident)))
library(gplots)
heatmap.2(HEATDATA,,scale=c("none"),dendrogram='none',Colv=T,trace='none',col=colorRampPalette(c('grey95','indianred')) ,margins=c(10,10))
library(dplyr)

used=which(!is.na(SUR[,5]))
boxplot(SUR[used,3],pch=16,ylab='OS (month)')$out
SUR[which(SUR[,3]==32),1]
SUR[which(SUR[,3]==82),1]


used=which((!is.na(SUR[,5])) & (!SUR[,3] %in% c(32,82) ))
SUR=read.table('SUR.txt',header=T)
boxplot(SUR[used,3]~SUR[used,5],pch=16,ylab='OS (month)')




library(survival)
library(survminer)

TYPE=SUR[used,5]
this_sur=SUR[used,3]

tmp=which(TYPE %in% c(1,2,3,4))
surtime=this_sur[tmp]
surevent=rep(1,length(tmp))
surtype=TYPE[tmp]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
ggsurvplot(fit, pval = TRUE)
surv_pvalue(fit)

used1=which(SUR[used,5]==1)
used2=which(SUR[used,5]==2)
used3=which(SUR[used,5]==3)
used4=which(SUR[used,5]==4)
t.test(SUR[used,2][used2],SUR[used,2][used4])


###########


pbmc@meta.data$newC= pbmc@meta.data$C
pbmc@meta.data$newC[which(names(pbmc@ident)=='Normal_PBTR.0010')]=10
pbmc@meta.data$newC[which(names(pbmc@ident)=='Tumor_PBTR.0010')]=10
pbmc@meta.data$newC[which(names(pbmc@ident)=='Normal_PBTR.0050')]=50
pbmc@meta.data$newC[which(names(pbmc@ident)=='Tumor_PBTR.0050')]=50


tmp=pbmc@ident
pbmc@ident=as.factor(pbmc@meta.data$newC)
names(pbmc@ident)=names(tmp)

#cluster2.markers <- FindMarkers(object = pbmc, ident.1 = 2, thresh.use = 0.25, 
#    test.use = "wilcox", only.pos = TRUE)
#write.table(cluster2.markers,file='C2_marker.txt',sep='\t',quote=F,row.names=T,col.names=T)
#pbmc.markers <- FindAllMarkers(object = pbmc,test.use = "bimod", only.pos = TRUE, min.pct = 0, thresh.use = 0.1)

#write.table(pbmc.markers,file='ALL_marker.txt',sep='\t',quote=F,row.names=T,col.names=T)
#library(dplyr)
#top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

#DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
this_used=which(!pbmc@meta.data$newC %in% c(10,50))
boxplot(pbmc@meta.data$nMut[this_used]~ pbmc@meta.data$newC[this_used])

c2_used=which(pbmc@meta.data$newC==2)

#c2v_mat=pbmc@raw.data[which(rownames(pbmc@raw.data) %in% pbmc@var.genes),c2_used]
c2v_mat=pbmc@raw.data[,c2_used]

c2v_mat_bin=as.matrix(c2v_mat)
c2v_mat_bin[which(c2v_mat_bin>0)]=1
c2v_mat_bin_num=apply(c2v_mat_bin,1,sum)


co_used=which(!pbmc@meta.data$newC %in% c(2,10,50))
cov_mat=pbmc@raw.data[,co_used]
cov_mat_bin=as.matrix(cov_mat)
cov_mat_bin[which(cov_mat_bin>0)]=1
cov_mat_bin_num=apply(cov_mat_bin,1,sum)


c2v_out=c2v_mat_bin_num[which(c2v_mat_bin_num>=2 & cov_mat_bin_num==0)]

which(rownames(c2v_mat)=='chr17.7578211.7578212.G>A')


write.table(c2v_out,file='C2V.txt',sep='\t',quote=F,row.names=T,col.names=F)
