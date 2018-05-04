
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("monocle")
#biocLite(c("DDRTree", "pheatmap"))
library(monocle)


MARKER=read.table('Marker_GENE.txt',header=T)
THIS=which(MARKER[,7] %in% c('Astro','RG','OPC') & MARKER[,6]<0.05)
MARKER_GENE=MARKER[THIS,1]




HSMM_expr_matrix <- read.table('2000_run1716_normalized-FOR QC-after qc_groupasseurat_only13689.txt.tsv',header=T,row.names=1)
HSMM_sample_sheet <- read.delim('2000_run1716_normalized-FOR QC-after qc_groupasseurat_only13689.txt.tsv.pheno',header=F,row.names=1)
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)

HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)


HSMM_myo <- setOrderingFilter(HSMM,  MARKER_GENE)
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)
plot_cell_trajectory(HSMM_myo, color_by = "Hours")


