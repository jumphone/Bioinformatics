library('Seurat')
GBM.data=read.table('./GBM/GBM_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
medullo.data=read.table('./medullo/medullo_out_gene_exon_tagged.dge.txt',header=T,row.names=1)
nb6.data=read.table('./nb6/nb6_out_gene_exon_tagged.dge.txt',header=T,row.names=1)


