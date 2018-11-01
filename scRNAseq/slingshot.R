source("https://bioconductor.org/biocLite.R")
biocLite("kstreet13/slingshot")

library(slingshot)
get_lineages(pcaX, clus, start.clus = 'm10')




#Author: Xinran Dong & Yuhao Feng
#refer to http://www.bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/slingshot.html#constructing-smooth-curves-and-ordering-cells
#refer to https://satijalab.org/seurat/pbmc3k_tutorial.html
library(Seurat)
library(dplyr)
library(slingshot)
library(RColorBrewer)
pbmc.data=read.table("/home/fyh/Desktop/slingshot/run1978_lane12_normalized_P_only OPC_tscan_marker_delete1290.txt",sep = "\t",header = TRUE,row.names = 1)
pbmc<-CreateSeuratObject(raw.data = pbmc.data,project = "10X_PBMC")
pbmc<-NormalizeData(object = pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:4, resolution = 0.2, print.output = 0, save.SNN = TRUE,force.recalc=TRUE)
t1   <- PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 6)
t2   <- cbind(t1$data$x,t1$data$y) 
rownames(t2) <- pbmc@cell.names 
g1=getLineages(t2,t1$data$ident,start.clus = 1)
#plot(t2,pch=16,col = brewer.pal(3,"Set1")[t1$data$ident],xlab='PC1',ylab='PC2')
#lines(g1)
crv1<-getCurves(g1)
crv1
plot(t2,pch=16,col = brewer.pal(3,"Set1")[t1$data$ident],xlab='PC1',ylab='PC2')
lines(crv1)
pseudo_order<-data.frame(slingPseudotime(crv1))
crv1_preorder<-data.frame(cbind(row.names(pseudo_order),pseudo_order$curve1))
crv1_order<-crv1_preorder[order(crv1_preorder$X2),]
write.table(crv1_order,file = "/home/fyh/Desktop/slingshot/crv1_order.txt",sep = "\t",row.names = FALSE)
# same for crv2




