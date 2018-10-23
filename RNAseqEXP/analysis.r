setwd("C:\\Data\\FDU_Hospital\\Liguo\\Control tumor and Olig2cKO tumor")
library("DOSE")
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)

library(gplots);
file="gene_matrix.txt"
data<-as.matrix(read.table(file,sep="\t",head=T,row.names=1,fill=TRUE))
 dim(data)
#[1] 11668    12
 colnames(data)
# [1] "gene"              "locus"             "sample_1"         
# [4] "sample_2"          "status"            "value_1"          
# [7] "value_2"           "log2.fold_change." "test_stat"        
#[10] "p_value"           "q_value"           "significant"  

col <-rep("Grey",nrow(data))
names(col)<-rownames(data)
data_s<-data[data[,10]<0.01,]
 dim(data_s)
#[1] 1809   12
p<--log(as.numeric(data[,10]))/log(10)
col[match(row.names(data_s[as.numeric(data_s[,8])>0,]),names(col))]<- "Red"
col[match(row.names(data_s[as.numeric(data_s[,8])<0,]),names(col))]<- "Blue"
pdf(file="volcano.pdf",10,10)
plot(data[,8],p,col=col, pch=19,xlab="Log2(fold change)",ylab="-Log10(p value)",xlim=c(-10,10),ylim=c(0,5),bty = 'n',pos=0)
dev.off()

gsea_dot<-read.table(file="GSEA_dotplot.txt",sep="\t",head=T,row.names=1)
dim(gsea_dot)
#[1] 34 10

p<- ggplot(gsea_dot[-c(1,2),], aes(x = NES, y = fct_reorder(GS.br..follow.link.to.MSigDB, NES))) + 
               geom_point(aes(color = FDR.q.val,size = SIZE)) +
               coord_cartesian(xlim=c(-3,3))+               
        scale_colour_gradientn(limits=c(0, 0.5), colours=rainbow(6)) +
        theme_bw(base_size = 14) +
        ylab(NULL) +scale_fill_brewer(palette = 'Accent')+theme(panel.background=element_rect(fill='grey95'),panel.border = element_blank())+
        geom_vline(xintercept = 0,size=1)+
        scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3))+
        ggtitle("GSEA analysis")

pdf(file="GSEA_dotplot.pdf",10,8) 
p
dev.off()     
        

pathway_matrix<-as.matrix(read.table(file="genelist_exp.txt",head=T,sep="\t",row.names=1))
 dim(pathway_matrix)
[1] 18  2


pdf ("heatmap_pathaway_mean.pdf",10,20)
heatmap.2(pathway_matrix/rowMeans(pathway_matrix),col=colorpanel(256, low="Blue", high="Red", mid="White"), density.info="none",trace="none")
dev.off()

