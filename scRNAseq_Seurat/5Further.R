library(Seurat)
library(dplyr)
library(Matrix)
load("./images/Seurat_EXP_TSNE.Robj")
load("./images/Seurat_EXP_cluster.Robj")


####Stat num of each cluster
i=0
while(i<=30){
CELL= which(EXP_cluster@ident==i)
TAB=as.matrix(table(EXP@ident[CELL]))
write.table(TAB,paste0('./images/NUM/NUM_',as.character(i),'.tsv'),sep='\t',row.names=T,col.names=F,quote=F)
i=i+1}
TAB=as.matrix(table(EXP@ident))
write.table(TAB,paste0('./images/NUM/NUM_total.tsv'),sep='\t',row.names=T,col.names=F,quote=F)
################################

#EXP@var.genes
C_2a=c(2,9,14,17,19,20,22,23)
C_2b=c(0,3,8,10,13,16,18,28)
C_2c=c(1,4,5,6,7,21)

exp_scale = t(apply(EXP@data,1,scale,scale=T,center=T))

Compare = function(CA,CB){
A=c()
B=c()
Z=c()
P=c()
i=1
while(i <=length(exp_scale[,1])){
a=mean(exp_scale[i,which(EXP_cluster@ident %in% CA)])
b=mean(exp_scale[i,which(EXP_cluster@ident %in% CB)])
z=a-b
pvalue=min(pnorm(z),1-pnorm(z))*2
A=c(A,a)
B=c(B,b)
Z=c(Z,z)
P=c(P,pvalue)
print(i)
i=i+1
}
return(list(A=A,B=B,Z=Z,P=P))
}

C_2a=c(2,9,14,17,19,20,22,23)
C_2b=c(0,3,8,10,13,16,18,28)
C_2c=c(1,4,5,6,7,21)


CCC=C_2c
OUT_DIR='./images/COMPARE/C_2c/'

i=1
while(i <=length(CCC)){
j=i+1
while(j <=length(CCC)){
Result=Compare(CCC[i],CCC[j])
H=c('GENE',as.character(CCC[i]),as.character(CCC[j]),'Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=paste0(OUT_DIR,as.character(CCC[i]),'_',as.character(CCC[j]),'.tsv'),quote=F,col.names=T,row.names=F,sep='\t')
j=j+1
print(j)
print(i)
}
i=i+1

}

Compare_age = function(CA,CB){
A=c()
B=c()
Z=c()
P=c()
i=1
while(i <=length(exp_scale[,1])){
a=mean(exp_scale[i,CA])
b=mean(exp_scale[i,CB])
z=a-b
pvalue=min(pnorm(z),1-pnorm(z))*2
A=c(A,a)
B=c(B,b)
Z=c(Z,z)
P=c(P,pvalue)
print(i)
i=i+1
}
return(list(A=A,B=B,Z=Z,P=P))
}

C_2a=c(2,9,14,17,19,20,22,23)
C_2b=c(0,3,8,10,13,16,18,28)
C_2c=c(1,4,5,6,7,21)


CCC=C_2a
TAG='C_2a'

CCC=C_2b
TAG='C_2b'


CCC=C_2c
TAG='C_2c'

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/1.5MD_4MD.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MD') 
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MD')
Result=Compare_age(CA,CB)
H=c('GENE','S1.5MD','S4MD','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/1.5MP_4MP.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MP') 
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MP')
Result=Compare_age(CA,CB)
H=c('GENE','S1.5MP','S4MP','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/1.5MSN_4MSN.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MSN') 
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MSN')
Result=Compare_age(CA,CB)
H=c('GENE','S1.5MSN','S4MSN','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')


OUT_FILE='./images/COMPARE_C/11_C2a.tsv'
CA=which(EXP_cluster@ident %in% c(11) )
CB=which(EXP_cluster@ident %in% C_2a)
Result=Compare_age(CA,CB)
H=c('GENE','11','C2a','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')



OUT_FILE='./images/COMPARE_C/11_C2b.tsv'
CA=which(EXP_cluster@ident %in% c(11) )
CB=which(EXP_cluster@ident %in% C_2b)
Result=Compare_age(CA,CB)
H=c('GENE','11','C2b','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')


OUT_FILE='./images/COMPARE_C/2022_C2a.tsv'
CA=which(EXP_cluster@ident %in% c(20,22) )
CB=which(EXP_cluster@ident %in% C_2a)
Result=Compare_age(CA,CB)
H=c('GENE','20_22','C2a','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')


###############################################################################################

TAG='C_2a'
CCC=C_2a



exp_scale = t(apply(EXP@data,1,scale,scale=T,center=T))


Compare_age = function(CA,CB){
A=c()
B=c()
Z=c()
P=c()
i=1
while(i <=length(exp_scale[,1])){
a=mean(exp_scale[i,CA])
b=mean(exp_scale[i,CB])
z=a-b
pvalue=min(pnorm(z),1-pnorm(z))*2
A=c(A,a)
B=c(B,b)
Z=c(Z,z)
P=c(P,pvalue)
print(i)
i=i+1}return(list(A=A,B=B,Z=Z,P=P))}





Compare_AGE=function(TAG,CCC){

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/1.5MD_1.5MP.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MD')
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MP')
Result=Compare_age(CA,CB)
H=c('GENE','S1.5MD','S1.5MP','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/1.5MD_1.5MSN.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MD')
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MSN')
Result=Compare_age(CA,CB)
H=c('GENE','S1.5MD','S1.5MSN','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/1.5MP_1.5MSN.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MP')
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S1.5MSN')
Result=Compare_age(CA,CB)
H=c('GENE','S1.5MP','S1.5MSN','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/4MD_4MP.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MD')
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MP')
Result=Compare_age(CA,CB)
H=c('GENE','S4MD','S4MP','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/4MD_4MSN.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MD')
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MSN')
Result=Compare_age(CA,CB)
H=c('GENE','S4MD','S4MSN','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')

OUT_FILE=paste0('./images/COMPARE_AGE/',TAG,'/4MP_4MSN.tsv')
CA=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MP')
CB=which(EXP_cluster@ident %in% CCC & EXP@ident %in% 'S4MSN')
Result=Compare_age(CA,CB)
H=c('GENE','S4MP','S4MSN','Z')
OUT=cbind(rownames(EXP_cluster@data),Result$A, Result$B, Result$Z)
colnames(OUT)=H
write.table(OUT,file=OUT_FILE,quote=F,col.names=T,row.names=F,sep='\t')



}


Compare_AGE('C_2a',C_2a)
Compare_AGE('C_2b',C_2b)
Compare_AGE('C_2c',C_2c)





i=0
while(i<=30){
print(i)
#TSNEPlot(object = EXP_cluster,colors.use=c(rep('grey',i),'red',rep('grey',30-i)))
cluster.markers <- FindMarkers(object = EXP_cluster, ident.1 = i, min.pct = 0.25)
cluser_top = head(x = cluster.markers, n = 1000)
write.table(file=paste0('./images/Cluster_1000/Cluster_',as.character(i),'_marker.tsv'),cluser_top,sep='\t',quote=F)
i=i+1}





#########################################################################################

library(Seurat)
library(dplyr)
library(Matrix)
load("./images/Seurat_EXP_TSNE.Robj")
load("./images/Seurat_EXP_cluster.Robj")
################################################################################################
C_2a=c(2,9,14,17,19,20,22,23)
C_2b=c(0,3,8,10,13,16,18,28)
C_2c=c(1,4,5,6,7,21)
library(TSCAN)
#EXP@var.genes

cell_cycle=read.table('CELL_CYCLE.txt.mouse')[,3]
TSCAN_data=EXP_cluster@data[which(!rownames(EXP_cluster@data) %in% cell_cycle) ,which(EXP_cluster@ident %in% C_2a )]
#TSCAN_data=EXP_cluster@data[,which(EXP_cluster@ident %in% C_2a )]
TSCAN_data=as.matrix(TSCAN_data)
rownames(TSCAN_data)=rownames(EXP_cluster@data)[which(!rownames(EXP_cluster@data) %in% cell_cycle)]
#rownames(TSCAN_data)=rownames(EXP_cluster@data)
colnames(TSCAN_data)=colnames(EXP_cluster@data)[which(EXP_cluster@ident %in% C_2a )]

VAR_ROW=apply(TSCAN_data,1,var)
AVE_ROW=apply(TSCAN_data,1,mean)
VAR_LIMIT=0.3
AVE_LIMIT=0.0
TSCAN_data_proc=TSCAN_data[which(VAR_ROW>VAR_LIMIT & AVE_ROW>AVE_LIMIT),]
rownames(TSCAN_data_proc)=rownames(TSCAN_data)[which(VAR_ROW>VAR_LIMIT & AVE_ROW>AVE_LIMIT)]

TSCAN_data_proc_clust=exprmclust(TSCAN_data_proc,clusternum=4, reduce = T,  cluster = NULL)
TSCAN_data_proc_clust_order=TSCANorder(TSCAN_data_proc_clust)

STEM=c('Lpar1','Prom1','Vim','Sox2','Sox9','Gdnf','Nes')
TSCAN_data_proc_stem = apply(TSCAN_data_proc[which(rownames(TSCAN_data_proc) %in% STEM),],2,sum)

pdf('./images/TSCAN_Pseudotime_singlegene.pdf',width=14,height=14)
names(TSCAN_data_proc_stem)= colnames(TSCAN_data_proc)
#plotmclust(TSCAN_data_proc_clust, x = 1, y = 2, MSTorder = NULL, show_tree = T,show_cell_names = F, cell_name_size = 3, markerexpr = TSCAN_data_proc_stem)
singlegeneplot(geneexpr=TSCAN_data_proc_stem,TSCAN_data_proc_clust_order,cell_size=2)
dev.off()


library(IDPmisc)
COLMP=colorRampPalette(c("grey98", "red"))
pdf('./images/TSCAN_PLOT.pdf',width=7,height=7)
plot(TSCAN_data_proc_clust$pcareduceres[,1],TSCAN_data_proc_stem,xlab='PC1',ylab='Stem_Gene_EXP',main=paste0('PCC=',as.character(cor(TSCAN_data_proc_clust$pcareduceres[,1],TSCAN_data_proc_stem))),pch=20)
abline(lm(TSCAN_data_proc_stem~TSCAN_data_proc_clust$pcareduceres[,1]),col='red')
#iplot(TSCAN_data_proc_clust$pcareduceres[,1],TSCAN_data_proc_stem,xlab='PC1',ylab='Stem_Gene_EXP',pixs=5,colramp=COLMP)
plot(TSCAN_data_proc_clust$pcareduceres[,2],TSCAN_data_proc_stem,xlab='PC2',ylab='Stem_Gene_EXP',main=paste0('PCC=',as.character(cor(TSCAN_data_proc_clust$pcareduceres[,2],TSCAN_data_proc_stem))),pch=20)
abline(lm(TSCAN_data_proc_stem~TSCAN_data_proc_clust$pcareduceres[,2]),col='red')
#iplot(TSCAN_data_proc_clust$pcareduceres[,2],TSCAN_data_proc_stem,xlab='PC2',ylab='Stem_Gene_EXP',pixs=5,colramp=COLMP)
plot(TSCAN_data_proc_clust$pcareduceres[,3],TSCAN_data_proc_stem,xlab='PC3',ylab='Stem_Gene_EXP',main=paste0('PCC=',as.character(cor(TSCAN_data_proc_clust$pcareduceres[,3],TSCAN_data_proc_stem))),pch=20)
abline(lm(TSCAN_data_proc_stem~TSCAN_data_proc_clust$pcareduceres[,3]),col='red')
#iplot(TSCAN_data_proc_clust$pcareduceres[,3],TSCAN_data_proc_stem,xlab='PC3',ylab='Stem_Gene_EXP',pixs=5,colramp=COLMP)
dev.off()

pdf('./images/TSCAN.pdf',width=30,height=30)
plotmclust(TSCAN_data_proc_clust, x = 1, y = 2, MSTorder = NULL, show_tree = T,show_cell_names = F, cell_name_size = 3, markerexpr = NULL)
plotmclust(TSCAN_data_proc_clust, x = 1, y = 3, MSTorder = NULL, show_tree = T,show_cell_names = F, cell_name_size = 3, markerexpr = NULL)
plotmclust(TSCAN_data_proc_clust, x = 2, y = 3, MSTorder = NULL, show_tree = T,show_cell_names = F, cell_name_size = 3, markerexpr = NULL)
dev.off()

##################################################################################################

OUT_FILE='./images/TSCAN/matrix_C2a.tsv'
write.table(TSCAN_data,file=OUT_FILE,quote=F,col.names=T,row.names=T,sep='\t')


pdf('images/Natalie_GENES_EXP.pdf',width=30,height=30)
FeaturePlot(object = EXP_cluster, features.plot = c("Plp1","Mbp",'Cnp','Pmp22','Sox10','Egr2','Pou3f1'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Lpar1","Prom1",'Vim','Sox2','Sox9','Gdnf','Nes'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Foxd3",'Tfap2','Ets1','Twist'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Cyr61",'Ankrd1','Amotl2','Ccnd1'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Pdgfa",'Pdgfb','Pdgfra','Pdgfrb'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Birc5",'Dusp6','Aldh1a1','Aurka','Top2a','Mki67'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Fbn1",'Fbn2','Spock1','Serpine1','Ets1','Fstl3','Fstl1','Edil3','Slit2','Plod2'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Csf1r",'Cxcl2','Cd68','Serpinh1'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Thy1",'Cd34','Fkbp7'), cols.use = c("grey", "blue"), reduction.use = "tsne")
FeaturePlot(object = EXP_cluster, features.plot = c("Cd63",'Cd81','Lamp1'), cols.use = c("grey", "blue"), reduction.use = "tsne")
dev.off()
