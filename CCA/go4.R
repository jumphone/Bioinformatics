

source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/BAST.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

D1=readRDS('CT.RDS')
D2=readRDS('MS.RDS')

D1=as.matrix(D1)
D2=as.matrix(D2)


bastout=BAST(D1,D2,FDR=0.05,COR=0.8)

DimPlot(object =bastout$seurat, reduction.use = "umap", group.by = "condition",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =bastout$seurat, reduction.use = "umap", group.by = "map",  pt.size = 0.1, do.return = TRUE)



library(scran)

gene.counts1=D1
sce1 <- SingleCellExperiment(list(counts=gene.counts1))
sce1 <- normalize(sce1)

gene.counts2=D2
sce2 <- SingleCellExperiment(list(counts=gene.counts2))
sce2 <- normalize(sce2)

b1 <- sce1
b2 <- sce2
out <- fastMNN(b1, b2)



