library(Seurat)
library(dplyr)
library(Matrix)


exp_data=read.table('data.txt.rmdup',header=T,row.names=1)

EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes=0)

allgene=rownames(EXP@data)

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

#EXP = ScaleData(object = EXP, vars.to.regress = c("nUMI"), genes.use = allgene)

EXP = ScaleData(object = EXP,  genes.use = allgene)

EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )

pdf('PCA.pdf',width=15,height=15)
PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)
dev.off()
