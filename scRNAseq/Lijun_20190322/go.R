library(Seurat)
pbmc.data <- read.table("GSE115011_marton_all_cells.csv.pure", sep=',', header=T,row.names=1)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200,  project = "10X_PBMC")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes) # 5126

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

PCNUM=100
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

PCElbowPlot(object = pbmc,num.pc=PCNUM)

PCUSE=1:100
pbmc <- RunUMAP(object = pbmc, dims.use = PCUSE, seed.use=1)

DimPlot(pbmc,reduction.use='umap')




#cell_use=colnames(pbmc@data)[which(pbmc@dr$umap@cell.embeddings[,2] > -5)]
#DimPlot(pbmc,reduction.use='umap',cells.use=cell_use)
pdf('f1_UMAP.pdf',width=7,height=8)
GENE=c("PDGFRA","CSPG4", "TOP2A",'MKI67','PCDH15' ,'MOG','MBP','CNP')
FeaturePlot(object = pbmc, features.plot = GENE, cols.use = c("grey", "blue"), reduction.use = "umap")
dev.off()






UMAP2=pbmc@dr$umap@cell.embeddings[,2]

V=(UMAP2-min(UMAP2))
V=V/max(V)
V=round(V*100)
V=V+1

COLP=colorRampPalette(c('blue3','darkgreen','gold','red3'))(101)

COL=COLP[V]


pdf('f3_TRA.pdf',width=5,height=10)
par(mfrow=c(5,1))
PDGFRA=pbmc@data[which(rownames(pbmc@data) == "PDGFRA"),]
plot(UMAP2, PDGFRA, pch=16,xlab='UMAP2',col=COL)
CSPG4=pbmc@data[which(rownames(pbmc@data) == "CSPG4"),]
plot(UMAP2, CSPG4, pch=16,xlab='UMAP2',col=COL)
TOP2A=pbmc@data[which(rownames(pbmc@data) == "TOP2A"),]
plot(UMAP2, TOP2A, pch=16,xlab='UMAP2',col=COL)
MBP=pbmc@data[which(rownames(pbmc@data) == "MBP"),]
plot(UMAP2, MBP, pch=16,xlab='UMAP2',col=COL)
CNP=pbmc@data[which(rownames(pbmc@data) == "CNP"),]
plot(UMAP2,CNP, pch=16,xlab='CNP',col=COL,ylim=c(0,10))
dev.off()


GS=c("ASCL1","PPP1R14B", "PDGFRA",'CSPG4','ZNF488' ,'NKX2-2','CNP','MYRF')
pdf('f5_TRA.pdf',width=5,height=2*length(GS))
par(mfrow=c(length(GS),1))
for(one in GS){
GENE=pbmc@data[which(rownames(pbmc@data) == one),]
plot(UMAP2, GENE, ylab=one, pch=16,xlab='UMAP2',col=COL)
}
dev.off()



pdf('f3_UMAP.pdf',width=7,height=8)
GENE=c("ASCL1","PPP1R14B", "PDGFRA",'CSPG4','ZNF488' ,'NKX2-2','CNP','MYRF')
FeaturePlot(object = pbmc, features.plot = GENE, cols.use = c("grey", "blue"), reduction.use = "umap")
dev.off()





saveRDS(pbmc,file='pbmc.RDS')





