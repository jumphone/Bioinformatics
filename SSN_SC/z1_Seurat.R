library(Seurat)
library(dplyr)
library(Matrix)

PCNUM=40
PCUSE=1:35
RES=0.8



load('../images/Seurat_EXP_cluster.Robj')

write.table(as.matrix(EXP_cluster@data),file='all_gene_data.txt',quote=F,sep='\t',row.names=T,col.names=T)



STDVAR=function(a){stdvar=sqrt(sum((a-mean(a))^2/(length(a)-1)));return(stdvar)}

exp_data=read.table('var_gene_data.txt.jaspar.result.tmp',header=T,row.names=1)

library(mice)
miceMod <- mice(exp_data, method="rf") 
miceOutput <- complete(miceMod)


#FACTOR=c()
#CNUM=length(exp_data[1,])
#TNUM=length(exp_data[,1])
#i=1
#while(i<=TNUM ){
#NANUM=length(which(is.na(t(exp_data[i,]))))
#NAR=NANUM/CNUM
#factor=(1-NAR)*(1-NAR)
#FACTOR=c(FACTOR, factor)
#i=i+1
#print(i)}


load("SSN_EXP.Robj")
FAC=read.table('tmp.result.factor',header=F)
new_data=EXP@data * FAC[,2]
EXP@scale.data=new_data
all_gene=row.names(EXP@data)
PCNUM=20
PCUSE=1:20
RES=0.3

EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE, check_duplicates = FALSE)
EXP_cluster <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE)

pdf('TF_SSN_Seurat.pdf',width=30,height=15)
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




pdf('PCA_A.pdf')
PCElbowPlot(object = EXP, num.pc = 20)
PrintPCA(object = EXP, pcs.print = 1:20)
PCHeatmap(object = EXP, pc.use = 1:10,100)
PCHeatmap(object = EXP, pc.use = 11:20,100)
dev.off()


new_exp_data = -log(exp_data,10)

exp_data[is.na(exp_data)]=0

EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=0)



stdvar=apply(exp_data,1,STDVAR)

all_gene=row.names(exp_data[which(stdvar>0.8),])
all_gene=row.names(exp_data)
                   
EXP = ScaleData(object = EXP,vars.to.regress = c("nUMI"), genes.use = all_gene)

EXP = ScaleData(object = EXP, genes.use = all_gene)

PCNUM=10
PCUSE=1:10
RES=0.6

EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE, check_duplicates = FALSE)
EXP_cluster <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE)

save(EXP, file = "SSN_EXP.Robj")
save(EXP_cluster, file = "SSN_EXP_cluster.Robj")

pdf('TF_SSN_Seurat.pdf',width=30,height=15)
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





