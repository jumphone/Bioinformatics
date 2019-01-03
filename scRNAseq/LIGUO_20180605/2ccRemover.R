#https://cran.r-project.org/web/packages/ccRemover/vignettes/ccRemover_tutorial.html

library(Seurat)
load('EXP.Robj')

#### ccRemover #############
library(Seurat)
library(ccRemover)
var_exp_data=as.matrix(EXP@data[which(rownames(EXP@data) %in% EXP@var.genes),])
mean_gene_exp <- rowMeans(var_exp_data)
exp_data_cen <- var_exp_data - mean_gene_exp 
gene_names <- rownames(exp_data_cen)
cell_cycle_gene_indices <- gene_indexer(gene_names, species = "mouse", name_type = "symbols" )
length(cell_cycle_gene_indices)
#1377
if_cc <- rep(FALSE,nrow(exp_data_cen)) 
if_cc[cell_cycle_gene_indices] <- TRUE
summary(if_cc)
exp_data_cen=as.matrix(exp_data_cen)
dat <- list(x=exp_data_cen, if_cc=if_cc)
xhat <- ccRemover(dat, bar=FALSE)
xhat <- xhat + mean_gene_exp
EXP@data[which(rownames(EXP@data) %in% EXP@var.genes),]=xhat
saveRDS(xhat,file='xhat.RDS')
############################


xhat=readRDS('xhat.RDS')

deNeg<-function(X){
   return(X-min(X))
}

pxhat=t(apply(xhat,1,deNeg))
colnames(pxhat)=colnames(xhat)
rownames(pxhat)=rownames(xhat)
pbmc = CreateSeuratObject(raw.data = pxhat, min.cells = 0, min.genes=0)
pbmc = NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc = ScaleData(object = pbmc)

EXP@scale.data=pbmc@scale.data


################################################################################################################




PCNUM=40
EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )


PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = EXP, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = EXP,num.pc=PCNUM)

PrintPCA(object = EXP, pcs.print = 1:20)



PCUSE=1:6
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE, check_duplicates = FALSE)

RES=0.5
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)





save(EXP, file = "EXP_ccRemover.Robj")

























