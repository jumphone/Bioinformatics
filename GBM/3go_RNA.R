#/home/zhangfeng/disk/project/ATAC/data


TAG=as.character(read.table('combined_RNA_rpkm_GBMTAG.txt',header=F)[,1])
a=read.table('combined_RNA_rpkm.txt',row.names=1,header=T,sep='\t',check.names = F)


library(Seurat)
library(dplyr)

pbmc <- CreateSeuratObject(raw.data = a, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc@meta.data$tag=TAG
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, do.plot=F,
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc)


PCNUM=15
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)

PCElbowPlot(object = pbmc)

PCUSE=1:15
pbmc <- RunTSNE(object = pbmc, perplexity=5, dims.use = PCUSE, do.fast = TRUE)

TSNEPlot(object = pbmc,group.by='tag',pt.size=4)
pdf('PCA.pdf')
PCAPlot(object = pbmc,group.by='tag',pt.size=4)
dev.off()
save.image(file='F1.RData')

################################

SURTAG=read.table('SURTAG.txt',sep='\t',header=F)

PC1=pbmc@dr$pca@cell.embeddings[,1]
PC1_GENE=pbmc@dr$pca@gene.loadings[,1]
GBM=which(TAG=='GBMx')
LGG=which(TAG=='LGGx')


D=which(SURTAG[,1]==1)
TAG_D=TAG[D]
COL=rep('red',length(D))
COL[which(TAG_D=='GBM')]='indianred3'
COL[which(TAG_D=='LGG')]='green3'

used=which(!is.na(SURTAG[,1]))
TAG_U=TAG[used]
COL_U=rep('red',length(used))
COL_U[which(TAG_U=='GBM')]='indianred3'
COL_U[which(TAG_U=='LGG')]='green3'

pdf('figure1.pdf')
plot(PC1[D],SURTAG[D,2],pch=16,col=COL, xlab='PC1',ylab='OS (days)', cex=1.5)
plot(PC1[used],SURTAG[used,2],pch=16,col=COL_U, xlab='PC1',ylab='OS (days)', cex=1.5)
dev.off()

cor.test(PC1[D],SURTAG[D,2],method='spearman')



library(survival)
library(survminer)
a=read.table('SUR.txt',header=T,sep='\t')
used=used
score=PC1[used]
SUR=SURTAG[used,2]
SURE=SURTAG[used,1]

pdf('figure2.pdf')
TYPE=rep('MED',length(used))
TYPE[which(score> quantile(score,0.5) )]='High'
TYPE[which(score<= quantile(score,0.5) )]='Low'
surtime=SUR[which(TYPE!='MED')]
surevent=SURE[(which(TYPE!='MED'))]
surtype=TYPE[which(TYPE!='MED')]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
gp=ggsurvplot(fit, pval = TRUE)
print(gp)
surv_pvalue(fit)
dev.off()

pdf('figure3.pdf')
plot(PC1[used],SURTAG[used,2],pch=16,col=COL_U, xlab='PC1',ylab='OS (days)', cex=1.5)
abline(v=quantile(score,0.5),col='red',lwd=1.5)
dev.off()
#############


#PrintPCA(object = pbmc, pcs.print = 1:1, genes.print = 10, use.full = FALSE)
PC1_PEAK_LOADING=pbmc@dr$pca@gene.loadings[,1]
PC2_PEAK_LOADING=pbmc@dr$pca@gene.loadings[,2]
PC3_PEAK_LOADING=pbmc@dr$pca@gene.loadings[,3]
PC4_PEAK_LOADING=pbmc@dr$pca@gene.loadings[,4]
HC= 0.005
LC= -0.005


par(mfrow=c(2,2))

plot(density(PC1_PEAK_LOADING))
abline(v=HC,col='red',lwd=1.5)
abline(v=LC,col='blue',lwd=1.5)

plot(density(PC2_PEAK_LOADING))
plot(density(PC3_PEAK_LOADING))
plot(density(PC4_PEAK_LOADING))

PC1_GENE=PC1_PEAK_LOADING
PC1_POS=PC1_GENE[which(PC1_GENE >HC)]
PC1_NEG=PC1_GENE[which(PC1_GENE <LC)]
length(PC1_POS)
length(PC1_NEG)

PC1_POS_N=names(PC1_POS)

tmp=strsplit(PC1_POS_N, "_")
PC1_POS_N_S=c()
for(one in tmp){
    PC1_POS_N_S = cbind(PC1_POS_N_S, one)
           }
PC1_POS_OUT=cbind(t(PC1_POS_N_S),PC1_POS)
colnames(PC1_POS_OUT)=c('CHR','START','END','TYPE','ID','PEAK_SCORE','REGION','GC','PC1_LOADING')
write.table(PC1_POS_OUT,'PC1_POS_OUT.txt',row.names=F,col.names=T,sep='\t',quote=F)

PC1_NEG_N=names(PC1_NEG)
tmp=strsplit(PC1_NEG_N, "_")
PC1_NEG_N_S=c()
for(one in tmp){
    PC1_NEG_N_S = cbind(PC1_NEG_N_S, one)
           }
PC1_NEG_OUT=cbind(t(PC1_NEG_N_S),PC1_NEG)

colnames(PC1_NEG_OUT)=c('CHR','START','END','TYPE','ID','PEAK_SCORE','REGION','GC','PC1_LOADING')
write.table(PC1_NEG_OUT,'PC1_NEG_OUT.txt',row.names=F,col.names=T,sep='\t',quote=F)


################################################################################################

PC1_POS_OUT=read.delim('PC1_POS_OUT.txt',sep='\t',header=T)
PC1_NEG_OUT=read.delim('PC1_NEG_OUT.txt',sep='\t',header=T)


rawbed = t(cbind(t(PC1_POS_OUT),t(PC1_NEG_OUT))) 
pc1_loading_scale= rep(2,nrow(rawbed))
pc1_loading_scale[which(as.numeric(rawbed[,9])<0)]= -2
#scale(as.numeric(rawbed[,9]))[,1]
bed=rawbed[,1:4]
colnames(bed)=c('chr','start','end','value')
bed=as.data.frame(bed)
bed$chr=as.character(rawbed[,1])
bed$start=as.numeric(rawbed[,2])
bed$end=as.numeric(rawbed[,3])
bed$value=pc1_loading_scale
bed=bed[order(bed[,2]),]
bed=bed[order(bed[,1]),]

UPLIMIT=2
LWLIMIT=-2
CEX=0.3
BCOL='grey30'
H=2

library(circlize)
circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", c(1:22, "X", "Y")))
circos.genomicTrackPlotRegion(bed, ylim = c(LWLIMIT, UPLIMIT), panel.fun = function(region, value, ...) {
    #print(value[,2])
   h=H
   cell.xlim = get.cell.meta.data("cell.xlim")
   #circos.lines(cell.xlim, c(h, h), col = BCOL)
        
   COL = rep('grey50',length(value[,1]))
   COL[which(value[,1]>0)]='red'
   COL[which(value[,1]<0)]='blue'
   #circos.genomicPoints(region, type='h',col=COL, value[,1], cex = CEX, pch = 16)
   X=region[,1]
   Y=value[,1]
   print(length(X))
   print(length(Y))
   circos.lines(x=X, type='h',col=COL, y=Y,baseline=0)
   
   #print(region)
   #circos.rect(region[1], 0, region[2], 0.5, col = 'red')
}, track.height = 0.2)



table(PC1_POS_OUT[,7])
table(PC1_NEG_OUT[,7])


par(mfrow=c(1,2))
pie(table(PC1_POS_OUT[,7]))
pie(table(PC1_NEG_OUT[,7]))


