library(Seurat)
library(dplyr)
library(Matrix)

load('SSN_EXP_cluster.Robj')
load('SSN_EXP.Robj')
SSN_CLUSTER=EXP_cluster@ident
SSN_CELL=colnames(EXP_cluster@data)

pdf('SSN_Seurat_Compare.pdf',width=20,height=10)
TSNEPlot(object = EXP_cluster,do.label=T)
plot1=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MD','S1.5MP','S1.5MSN','S4MD','S4MP','S4MSN'))
plot2=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MP','S1.5MD','S1.5MSN','S4MD','S4MP','S4MSN'))
plot3=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MSN','S1.5MD','S1.5MP','S4MD','S4MP','S4MSN'))
plot4=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MD','S1.5MP','S1.5MSN','S1.5MD','S4MP','S4MSN'))
plot5=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MP','S1.5MD','S1.5MSN','S4MD','S1.5MP','S4MSN'))
plot6=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MSN','S1.5MD','S1.5MP','S4MD','S4MP','S1.5MSN'))
plot_grid(plot1, plot2,plot3,plot4,plot5,plot6)


load('../../images/Seurat_EXP_cluster.Robj')
load('../../images/Seurat_EXP_TSNE.Robj')
TSNEPlot(object = EXP_cluster,do.label=T)
plot1=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MD','S1.5MP','S1.5MSN','S4MD','S4MP','S4MSN'))
plot2=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MP','S1.5MD','S1.5MSN','S4MD','S4MP','S4MSN'))
plot3=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MSN','S1.5MD','S1.5MP','S4MD','S4MP','S4MSN'))
plot4=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MD','S1.5MP','S1.5MSN','S1.5MD','S4MP','S4MSN'))
plot5=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MP','S1.5MD','S1.5MSN','S4MD','S1.5MP','S4MSN'))
plot6=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MSN','S1.5MD','S1.5MP','S4MD','S4MP','S1.5MSN'))
plot_grid(plot1, plot2,plot3,plot4,plot5,plot6)


OLD_CLUSTER=EXP_cluster@ident
OLD_CELL=colnames(EXP_cluster@data)
TMP=as.character(OLD_CLUSTER)
SSN_CLUSTER=as.character(SSN_CLUSTER)
TMP[which(OLD_CELL %in% SSN_CELL)]=SSN_CLUSTER
NAF=as.character(max(as.numeric(SSN_CLUSTER))+1)
TMP[which(! OLD_CELL %in% SSN_CELL)]=NAF
TMP=as.factor(TMP)
names(TMP)=names(EXP_cluster@ident)
EXP_cluster@ident=TMP


TSNEPlot(object = EXP_cluster,do.label=T)
plot1=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MD','S1.5MP','S1.5MSN','S4MD','S4MP','S4MSN'))
plot2=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MP','S1.5MD','S1.5MSN','S4MD','S4MP','S4MSN'))
plot3=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S1.5MSN','S1.5MD','S1.5MP','S4MD','S4MP','S4MSN'))
plot4=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MD','S1.5MP','S1.5MSN','S1.5MD','S4MP','S4MSN'))
plot5=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MP','S1.5MD','S1.5MSN','S4MD','S1.5MP','S4MSN'))
plot6=TSNEPlot(object = EXP,do.return = T,colors.use=c('grey','grey','grey','grey','grey','red'),plot.order=c('S4MSN','S1.5MD','S1.5MP','S4MD','S4MP','S1.5MSN'))
plot_grid(plot1, plot2,plot3,plot4,plot5,plot6)


ALL=c(0:(length(table(EXP_cluster@ident))-1))

i=0
while(i<length(table(EXP_cluster@ident))){
print(i)
TSNEPlot(object = EXP_cluster,do.return = F,colors.use=c(rep('grey',length(table(EXP_cluster@ident))-1),'red'), plot.order=as.character(c(i,ALL[which(! ALL %in% i)])) )
i=i+1
}

ALL=c(1:(length(table(EXP_cluster@ident))))

i=1
while(i<=length(table(EXP_cluster@ident))){
print(i)
TSNEPlot(object = EXP_cluster,do.return = F,colors.use=c(rep('grey',length(table(EXP_cluster@ident))-1),'red'), plot.order=as.character(c(i,ALL[which(! ALL %in% i)])) )
i=i+1
}



dev.off()




load('../images/Seurat_EXP_cluster.Robj')
#EXP_cluster@data=EXP@data
pdf('TF_SIG_ORI.pdf',width=20,height=10)
i=1
while(i<=length(EXP@data[,1])){
tf_name=rownames(EXP@data)[i]
FeaturePlot(object = EXP_cluster, features.plot = c(tf_name), cols.use = c("grey", "red"), reduction.use = "tsne")
i=i+1}
dev.off()


