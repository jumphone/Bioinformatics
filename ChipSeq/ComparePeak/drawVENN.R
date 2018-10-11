# install.packages('VennDiagram')
args<-commandArgs(T)
library(VennDiagram)


a=read.table(paste0(args[1],'/compare.txt'),sep='\t',header=F)
A=a[which(a[,5]>0),4]
B=a[which(a[,6]>0),4]

venn.diagram(x=list(F1=A, F2=B ), paste0(args[1],"/VENN.png"), height = 450, width = 450, resolution =300, imagetype="png", col="white", fill=c(colors()[616], colors()[38]), alpha=c(0.6, 0.6),lwd=0.5, cex=0.5,cat.cex=0.5)
