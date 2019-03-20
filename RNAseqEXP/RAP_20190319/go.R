library(Seurat)
library(dplyr)


pbmc.data <-read.table('combined.txt',sep='\t',row.names=1,header=T)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))


V=apply(pbmc@data,1,var)
    
used_gene=rownames(pbmc@data)[which(V>=quantile(V, 0.8))]
pbmc <- RunPCA(object = pbmc, pc.genes = used_gene, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc@meta.data$name=colnames(pbmc@data)



pdf('PCA08.pdf',width=14,height=12)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2,group.by='name', do.label=T)
dev.off()




    
used_gene=rownames(pbmc@data)[which(V>0)]



pbmc <- RunPCA(object = pbmc, pc.genes = used_gene, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc@meta.data$name=colnames(pbmc@data)



pdf('PCA_ALL.pdf',width=14,height=12)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2,group.by='name', do.label=T)
dev.off()
