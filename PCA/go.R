library(Seurat)
library(dplyr)
library(Matrix)
library(htmlwidgets)
library(plotly)


exp_data=read.table('data.txt.rmdup',header=T,row.names=1)
##########################################

EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes=0)

allgene=rownames(EXP@data)

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

#EXP = ScaleData(object = EXP, vars.to.regress = c("nUMI"), genes.use = allgene)

EXP = ScaleData(object = EXP,  genes.use = allgene)

PCNUM=3
EXP <- RunPCA(object = EXP, pc.genes = allgene, do.print = TRUE, pcs.print = 1:3,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )

OUT=PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2, do.label=F, no.legend =T, do.hover = T)
htmlwidgets::saveWidget(OUT, file = "index.html")

##########################################


exp_data=read.table('data.txt.rmdup',header=T,row.names=1,check.names =F)
exp_data[which(is.na(exp_data))]=0
exp_data_scale=apply(exp_data,1,scale)
exp_data_sum=apply(exp_data,1,sum)
exp_data_var=apply(exp_data,1,var)
exp_data_scale_t=t(exp_data_scale)
rownames(exp_data_scale_t)=rownames(exp_data)
colnames(exp_data_scale_t)=colnames(exp_data)
exp_data_scale_t=exp_data_scale_t[which(exp_data_var>0),]
exp_data_scale_t=as.matrix(exp_data_scale_t)
exp_data_scale_t[which(is.na(exp_data_scale_t))]=0
exp_data_scale_t_pca=princomp(exp_data_scale_t)
set.seed(123)

PC1=exp_data_scale_t_pca$loadings[,1]
PC2=exp_data_scale_t_pca$loadings[,2]
PC1=jitter(PC1,factor=50)
PC2=jitter(PC2,factor=50)
LABEL=rownames(exp_data_scale_t_pca$loadings)
pdf('PCA.pdf',width=7,height=7)
plot(PC1,PC2,pch=16, cex=3,col=c('green3','green3','red','red','black','black','blue','blue'))
text(PC1, PC2, labels=LABEL, cex= 0.7, pos=c(1,2,3,4))
dev.off()
