
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("monocle")
#biocLite(c("DDRTree", "pheatmap"))
library(monocle)


MARKER=read.table('Marker_GENE.txt',header=T)
THIS=which(MARKER[,7] %in% c('Astro','RG','OPC','IGC') & MARKER[,6]<0.05)
#THIS=which( MARKER[,6]<0.05)
MARKER_GENE=MARKER[THIS,1]
#MARKER_GENE=c('Olig2','Myc')



HSMM_expr_matrix <- read.table('2000_run1716_normalized-FOR QC-after qc_groupasseurat_only135689.txt.tsv',header=T,row.names=1)
HSMM_sample_sheet <- read.delim('2000_run1716_normalized-FOR QC-after qc_groupasseurat_only135689.txt.tsv.pheno',header=F,row.names=1)
colnames(HSMM_sample_sheet)=c('LABEL')
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)


HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

print(head(pData(HSMM)))

HSMM_myo <- setOrderingFilter(HSMM,  MARKER_GENE)
plot_ordering_genes(HSMM_myo)

#HSMM_myo <- reduceDimension(HSMM_myo, max_components = 3,method = 'DDRTree')
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 3,method = 'tSNE')
HSMM_myo <- orderCells(HSMM_myo)
plot_cell_trajectory(HSMM_myo)

HSMM_myo <- orderCells(HSMM_myo,root_state=4)
plot_cell_trajectory(HSMM_myo,color_by = "LABEL")


