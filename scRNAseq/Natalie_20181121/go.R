library('Seurat')

##########################################
GBM.data=read.table('./GBM/GBM_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
GBM <- CreateSeuratObject(raw.data = GBM.data, min.cells = 3, min.genes = 200, project = "GBM")

##########################################
medullo.data=read.table('./medullo/medullo_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
medullo <- CreateSeuratObject(raw.data = medullo.data, min.cells = 3, min.genes = 200, project = "medullo")

##########################################
nb6.data=read.table('./nb6/nb6_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
nb6 <- CreateSeuratObject(raw.data = nb6.data, min.cells = 3, min.genes = 200, project = "nb6")

