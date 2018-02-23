library(Seurat)
library(dplyr)
library(Matrix)


#S1.5MD=tag[which(tag[,2]=='S1.5MD'),1]
#S1.5MP=tag[which(tag[,2]=='S1.5MP'),1]
#S1.5MSN=tag[which(tag[,2]=='S1.5MSN'),1]
#S4MSN=tag[which(tag[,2]=='S4MSN'),1]
#S4MP=tag[which(tag[,2]=='S4MP'),1]
#S4MD=tag[which(tag[,2]=='S4MD'),1]


PCNUM=40
PCUSE=1:35
RES=2

exp_data=read.table('../data_new/EXP123456.combined.txt',header=T,row.names=1)
EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=200)

mito.genes <- grep(pattern = "^Mt", x = rownames(x = EXP@data), value = TRUE)
percent.mito <- colSums(EXP@data[mito.genes, ]) / colSums(EXP@data)
EXP <- AddMetaData(object = EXP, metadata = percent.mito, col.name = "percent.mito")


pdf('./mages/Seurat_QC.pdf',width=30,height=15)
VlnPlot(object = EXP, features.plot = c("nGene", "nUMI",'percent.mito'), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "nGene")
dev.off()

EXP=FilterCells(object = EXP, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(3000, 0.05))

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

pdf('./images/Seurat_VarGene.pdf')
EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.8)
dev.off()

length(x=EXP@var.genes)

#EXP = ScaleData(object = EXP,vars.to.regress = c("percent.mito", "orig.ident", "nUMI"), genes.use = EXP@var.genes, model.use = "negbinom")
EXP = ScaleData(object = EXP,vars.to.regress = c("percent.mito", "nUMI"), genes.use = EXP@var.genes)
EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )

PrintPCA(object = EXP, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

pdf('./images/Seurat_PCA.pdf')
PCElbowPlot(object = EXP,num.pc=PCNUM)
VizPCA(object = EXP, pcs.use = 1:5)
PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)
PCAPlot(object = EXP, dim.1 = 1, dim.2 = 3)
PCAPlot(object = EXP, dim.1 = 1, dim.2 = 4)
PCAPlot(object = EXP, dim.1 = 1, dim.2 = 5)
dev.off()

pdf('./images/PCA_A.pdf')

PCElbowPlot(object = EXP, num.pc = 40)
PrintPCA(object = EXP, pcs.print = 1:36)
PCHeatmap(object = EXP, pc.use = 1:12,100)
PCHeatmap(object = EXP, pc.use = 13:24,100)
PCHeatmap(object = EXP, pc.use = 25:36,100)
PCHeatmap(object = EXP, pc.use = 37:40,100)
dev.off()

save(EXP, file = "./images/Seurat_EXP_PCA.Robj")



EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE)
save(EXP, file = "./images/Seurat_EXP_TSNE.Robj")


EXP_cluster <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE)
save(EXP_cluster, file = "./images/Seurat_EXP_cluster.Robj")

pdf('./images/Seurat_TSNE.pdf',width=30,height=15)
TSNEPlot(object = EXP)
TSNEPlot(object = EXP_cluster)
TSNEPlot(object = EXP_cluster,do.label=T)
plot1=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MD','S1.5MP','S1.5MSN','S4MD','S4MP','S4MSN'))
plot2=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MP','S1.5MD','S1.5MSN','S4MD','S4MP','S4MSN'))
plot3=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MSN','S1.5MD','S1.5MP','S4MD','S4MP','S4MSN'))
plot4=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MD','S1.5MP','S1.5MSN','S1.5MD','S4MP','S4MSN'))
plot5=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MP','S1.5MD','S1.5MSN','S4MD','S1.5MP','S4MSN'))
plot6=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MSN','S1.5MD','S1.5MP','S4MD','S4MP','S1.5MSN'))
plot_grid(plot1, plot2,plot3,plot4,plot5,plot6)
dev.off()


pdf('images/GENES_EXP.pdf',width=30,height=30)
FeaturePlot(object = EXP_cluster, features.plot = c("Adar","Cntf",'Sox10','Gfap','Map2','Dlg3','Serpini1','Bmp2','S100b','Aldh1l1','Aif1','Alk'), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()


pdf('./images/Seurat_TSNE_clusters.pdf',width=30,height=15)
#TSNEPlot(object = EXP_cluster,colors.use=c(rep('grey',11),'red',rep('grey',19)))
i=0
while(i<=30){
print(i)
TSNEPlot(object = EXP_cluster,colors.use=c(rep('grey',i),'red',rep('grey',30-i)))
cluster.markers <- FindMarkers(object = EXP_cluster, ident.1 = i, min.pct = 0.25)
cluser_top = head(x = cluster.markers, n = 1000)
write.table(file=paste0('./images/Cluster_1000/Cluster_',as.character(i),'_marker.tsv'),cluser_top,sep='\t',quote=F)
i=i+1}
dev.off()



pdf('test.pdf',width=30,height=30)
FeaturePlot(object = EXP, features.plot = c('Pmp2','Mme','Fam107a'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP, features.plot = c('Ccnd1','Gstm1','Ankrd1'), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
