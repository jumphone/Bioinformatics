
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


exp_mat=read.delim('matt_combined.txt',sep='\t',row.names=1,header=T)

out=.get_cor(exp_mat, exp_mat, method='spearman',CPU=4)



library('gplots')

out[which(out<0)]=0
heatmap.2(out,scale=c("none"),dendrogram='both',Colv=T,trace='none',col=colorRampPalette(c('white','yellow','red1','red2','red3')) ,margins=c(15,15))







library(Seurat)
pbmc.data=read.delim('matt_combined.txt',sep='\t',row.names=1,header=T)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- ScaleData(object = pbmc,vars.to.regress = c("nUMI"))

PC1_POS=as.character(read.delim('PC1_POS_OUT.txt.com',header=F,sep='\t')[,1])
PC1_NEG=as.character(read.delim('PC1_NEG_OUT.txt.com',header=F,sep='\t')[,1])
PC1P=which(rownames(pbmc@scale.data) %in% PC1_POS)
PC1N=which(rownames(pbmc@scale.data) %in% PC1_NEG)

getPC1score=function(X){
    return(mean(X[PC1P])-mean(X[PC1N]))

    }

PC1score=apply(pbmc@scale.data,2,getPC1score)
par(mar=c(15,5,5,5))
barplot(PC1score,las=2,ylab='PC1 Score')


##########################
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI"), nCol = 2)

pbmc <- FilterCells(object = pbmc, subset.names = c("nUMI"), low.thresholds = c(1000000))




pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)

length(x = pbmc@var.genes)



PCNUM=10
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, pcs.compute=PCNUM, do.print = TRUE, pcs.print = 1:5, genes.print = 5)


PCUSE=1:10
pbmc <- RunTSNE(object = pbmc, perplexity=5, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE)


pbmc@meta.data$tag=colnames(pbmc@data)
TSNEPlot(object = pbmc,pt.size=4, do.label=T)
TSNEPlot(object = pbmc,pt.size=4,group.by='tag', do.label=T)
#PCAPlot(object = pbmc,pt.size=4, do.label=T)
