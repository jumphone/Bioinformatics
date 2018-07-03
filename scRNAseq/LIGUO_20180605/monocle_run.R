#source("https://bioconductor.org/biocLite.R")
#biocLite("monocle")
#biocLite(c("DDRTree", "pheatmap"))
#expectation   <- data.frame(expectation=model_expectation[cds_exprs$f_id[1], cds_exprs$Cell])

library(monocle)
library(plyr)
library(ggplot2)

exp_data=read.table('SUBDATA_var_500.txt',header=T,row.names=1)

exp_data=as.matrix(exp_data)
all_genes=rownames(exp_data) 
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
remove(HSMM)

ordering_genes = all_genes
#ordering_genes = c('Olig2','Olig1')

HSMM_myo <- estimateSizeFactors(HSMM_myo)
HSMM_myo = estimateDispersions(HSMM_myo)
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, reduction_method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)

plot_cell_trajectory(HSMM_myo, color_by = "State")
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")

show_genes <- row.names(subset(fData(HSMM_myo), gene_short_name %in% c("Olig1")))

plot_genes_jitter(HSMM_myo[show_genes,], grouping = "State", min_expr = 0.1, color_by = "State")




