
library(Seurat)
library(dplyr)
library(Matrix)
load('../../images/Seurat_EXP_cluster.Robj')

a=read.table('test0.2.data.result.b.KM_100.K20_TMP/tmp.20.k')
a=as.factor(t(a))
names(a)=names(EXP_cluster@ident)
EXP_cluster@ident=a
pdf('NEW.pdf',width=20,height=10)
TSNEPlot(object = EXP_cluster,do.label=T)
ALL=c(0:(length(table(EXP_cluster@ident))-1))
i=0
while(i<length(table(EXP_cluster@ident))){
print(i)
TSNEPlot(object = EXP_cluster,do.return = F,colors.use=c(rep('grey',length(table(EXP_cluster@ident))-1),'red'), plot.order=as.character(c(i,ALL[which(! ALL %in% i)])) )
i=i+1
}
dev.off()

