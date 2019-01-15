oexp=read.table('run1978lane12_out_gene_exon_tagged.dge_priopciol.txt.uniq.txt',sep='\t',header=F)
otag=as.character(t(oexp[1,])[2:ncol(oexp),])
odata=oexp[3:nrow(oexp),2:ncol(oexp)]
odata=apply(odata,2,as.numeric)
rownames(odata)=oexp[3:nrow(oexp),1]
colnames(odata)=as.character(t(oexp[2,])[2:ncol(oexp),])
odata=as.matrix(odata)
source('scRef.R')
our_data=.generate_ref(odata,cbind(otag,otag))


ref_tag=read.table('EXPMAT_TAG.txt',sep='\t',header=F)


rtag=as.character(ref_tag[,4])
rdata=read.table('EXPMAT.txt.uniq.txt',row.names=1,header=T,sep='\t')
ref_data=.generate_ref(rdata,cbind(rtag,rtag))


Pdgfra_our_index=which(rownames(odata)=='Pdgfra')
Pdgfra_our=odata[Pdgfra_our_index,]
Pdgfra_our[which(Pdgfra_our>0)]=1
table(Pdgfra_our,otag)
t(table(Pdgfra_our,otag))/apply(table(Pdgfra_our,otag),2,sum)

Pdgfra_ref_index=which(rownames(rdata)=='Pdgfra')
Pdgfra_ref=rdata[Pdgfra_ref_index,]
Pdgfra_ref[which(Pdgfra_ref>0)]=1
table(Pdgfra_ref,rtag)
t(table(Pdgfra_ref,rtag))/apply(table(Pdgfra_ref,rtag),2,sum)











log2p1=function(X){return(log(X+1,2))}
our_data_log=apply(our_data,2,log2p1)
ref_data_log=apply(ref_data,2,log2p1)
colnames(ref_data_log)=c('preOPC','OPC')

our_data_log=our_data_log
out=.get_cor(our_data_log,ref_data_log,method='pearson')
library('gplots')
heatmap.2(out,scale='col',margin=c(15,15),dendrogram='none',trace='none',col=colorRampPalette(c('blue','white','red')))


plot(ref_tag[,2],ref_tag[,3],col=as.factor(ref_tag[,4]),pch=16,xlab='tSNE1',ylab='tSNE2',
	xlim=c(-10,-2),ylim=c(-20,-12))

pdf('tSNEori.pdf',width=4.5,height=5)
plot(ref_tag[,2],ref_tag[,3],col=as.factor(ref_tag[,4]),pch=16,xlab='tSNE1',ylab='tSNE2',
	xlim=c(-10,-2),ylim=c(-20,-11))
dev.off()

exp_sc_mat=our_data
exp_ref_mat=ref_data
exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
gene_sc=rownames(exp_sc_mat)
gene_ref=rownames(exp_ref_mat)
gene_over= gene_sc[which(gene_sc %in% gene_ref)]
exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
colname_sc=colnames(exp_sc_mat)
colname_ref=colnames(exp_ref_mat)


colnames(exp_ref_mat)=c('bpreOPC','bOPC')
#AD=cbind(exp_sc_mat[,2:3],exp_ref_mat)
AD=cbind(exp_sc_mat,exp_ref_mat)
library(Seurat)
library(dplyr)

pbmc <- CreateSeuratObject(raw.data = AD, min.cells = 3, min.genes = 0, 
    project = "10X_PBMC")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)


pbmc <- ScaleData(object = pbmc, vars.to.regress = c('nUMI', "percent.mito"))
pbmc <- FindVariableGenes(object = pbmc,do.plot=F, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.18, x.high.cutoff = 2.25, y.cutoff = 0.5)
length(x = pbmc@var.genes)
allgene=rownames(pbmc@data)
pbmc@meta.data$tag=colnames(pbmc@data)
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, pcs.compute=3, do.print = F)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2,group.by='tag',pt.size=6)



log2p1=function(X){return(log(X+1,2))}
our_data_log= pbmc@scale.data[,1:3]#apply(our_data,2,log2p1)
ref_data_log= pbmc@scale.data[,4:5]#apply(ref_data,2,log2p1)
colnames(ref_data_log)=c('preOPC','OPC')
our_data_log=our_data_log
out=.get_cor(our_data_log,ref_data_log,method='pearson')


priOPC=pbmc@scale.data[,3]
preOPC=pbmc@scale.data[,4]
cor.test(priOPC,preOPC)

plot(priOPC,preOPC)


pdf('PCA_PC1.pdf',width=5,height=4)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2,group.by='tag',pt.size=6)
dev.off()




pbmc <- RunTSNE(object = pbmc, dims.use = 1:3, do.fast = F,perplexity=1)
TSNEPlot(object = pbmc,group.by='tag',pt.size=10)


pdf('COR.pdf',width=5,height=5)
plot(log(exp_sc_mat[,3]+1),log(exp_ref_mat[,1]+1),cex=0.2,
	pch=16,ylab=c("Gene's expression - preOPC"),xlab=c("Gene's expression - priOPC"))
dev.off()


priOPC=log(exp_sc_mat[,3]+1)
preOPC=log(exp_ref_mat[,1]+1)

names(which(priOPC>6 & preOPC<3))

cor.test(preOPC,priOPC)

PC1=pbmc@dr$pca@cell.embeddings[,1]
Pdgfra=pbmc@data[which(allgene=='Pdgfra'),]
Pax6=pbmc@data[which(allgene=='Pax6'),]
Gfap=pbmc@data[which(allgene=='Gfap'),]

par(mfrow=c(1,2))
plot(PC1, Pdgfra, pch=16,col=as.factor(names(PC1)))
plot(PC1, Pax6, pch=16,col=as.factor(names(PC1)))
#plot(PC1, Gfap, pch=16,col=as.factor(names(PC1)))


