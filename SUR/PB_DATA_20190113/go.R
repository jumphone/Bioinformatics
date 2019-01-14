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
TSNEPlot(object = pbmc,pt.size=3)




