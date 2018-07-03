#source("https://bioconductor.org/biocLite.R")
#biocLite("monocle")
#biocLite(c("DDRTree", "pheatmap"))

library(monocle)

exp_data=read.table('SUBDATA_var_500.txt',header=T,row.names=1)

exp_data=as.matrix(exp_data)
col_data=colnames(exp_data)
names(col_data)=colnames(exp_data)
col_data=as.data.frame(col_data)
row_data=rownames(exp_data)
names(row_data)=rownames(exp_data)
row_data=as.matrix(row_data)
colnames(row_data)='gene_short_name'
row_data=as.data.frame(row_data)

pd <- new("AnnotatedDataFrame", data =col_data)
fd <- new("AnnotatedDataFrame", data =row_data)
HSMM=newCellDataSet(exp_data, phenoData = pd, featureData = fd)
HSMM_myo=HSMM

diff_test_res <- differentialGeneTest(HSMM_myo)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
#ordering_genes = c('Olig2','Olig1')

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)

plot_cell_trajectory(HSMM_myo)
plot_cell_trajectory(HSMM_myo, color_by = "State")
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
