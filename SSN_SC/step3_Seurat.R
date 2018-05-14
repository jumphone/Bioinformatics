library(Seurat)
library(dplyr)
library(Matrix)

PCNUM=40
PCUSE=1:35
RES=0.8


exp_data=read.table('OUT',header=T,row.names=1)
exp_data[is.na(exp_data)]=0
exp_data=10^exp_data

EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=200)
EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)
pdf('./VarGene.pdf')
EXP=FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.8)
dev.off()
length(x=EXP@var.genes)

all_gene=EXP@var.genes #row.names(exp_data)
EXP = ScaleData(object = EXP,vars.to.regress = c("nUMI"), genes.use = all_gene)
EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE)
EXP_cluster <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE)

save(EXP, file = "NEW_SSN_EXP.Robj")
save(EXP_cluster, file = "NEW_SSN_EXP_cluster.Robj")

pdf('NEW_SSN_Seurat.pdf',width=30,height=15)
TSNEPlot(object = EXP)
TSNEPlot(object = EXP_cluster)
TSNEPlot(object = EXP_cluster,do.label=T)
plot1=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MD','S1.5MP','S1.5MSN','S4MD','S4MP','S4MSN'))
plot2=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MP','S1.5MD','S1.5MSN','S4MD','S4MP','S4MSN'))
plot3=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MSN','S1.5MD','S1.5MP','S4MD','S4MP','S4MSN'))
plot4=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MD','S1.5MP','S1.5MSN','S1.5MD','S4MP','S4MSN'))
plot5=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MP','S1.5MD','S1.5MSN','S4MD','S1.5MP','S4MSN'))
plot6=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MSN','S1.5MD','S1.5MP','S4MD','S4MP','S1.5MSN'))
plot_grid(plot1, plot2,plot3,plot4,plot5,plot6)
dev.off()





