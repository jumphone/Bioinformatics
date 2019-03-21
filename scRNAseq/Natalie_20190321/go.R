library(Seurat)
load('Lats12SN2samples.RData')
DimPlot(GFPSN_merged, reduction.use='tsne', group.by='orig.ident')




