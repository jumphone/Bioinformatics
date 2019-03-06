
source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/BEER.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

D1=readRDS('CT.RDS')
D2=readRDS('MS.RDS')


beerout=BEER(D1, D2, CNUM=100, PCNUM=50, CPU=4, print_step=10)


pbmc=beerout$seurat

pbmc <- RunUMAP(object = pbmc, reduction.use='adjpca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
DimPlot(pbmc,reduction.use='umap',group.by='condition',pt.size=0.1)
DimPlot(pbmc,reduction.use='umap',group.by='map',pt.size=0.1)

pcpbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
DimPlot(pcpbmc,reduction.use='umap',group.by='condition',pt.size=0.1)
DimPlot(pcpbmc,reduction.use='umap',group.by='map',pt.size=0.1)













DR=beerout$seurat@dr$pca@cell.embeddings
GROUP=c(beerout$g1,beerout$g2)
B1index=c(1:length(beerout$g1))
B2index=c((length(beerout$g1)+1):(length(beerout$g1)+length(beerout$g2)))
VP=beerout$vp
VPC=beerout$vpcor

.dr2adr <- function(DR, B1index, B2index, GROUP, VP){
    #set.seed(SEED)
    #library(dtw)
    #library(MALDIquant)
    #library(pcaPP)
    OUT=list()
    OUT$adr=DR
    VALID_PAIR=VP
    ALL_COR=c()   
    ALL_PV=c() 
    ALL_UCOR=c()   
    ALL_UPV=c() 
    index1=B1index
    index2=B2index
  
    vindex1=which(GROUP %in% VP[,1])
    vindex2=which(GROUP %in% VP[,2])
  
    
    print('Start')
    THIS_DR=1
    while(THIS_DR<=ncol(DR)){
        THIS_PC = DR[,THIS_DR]
       
        ########################
        sd_lst1=c()
        mean_lst1=c()
        max_lst1=c()
        min_lst1=c()
        sd_lst2=c()
        mean_lst2=c()  
        max_lst2=c()
        min_lst2=c()
        
        
        
        i=1
        while(i<=nrow(VP)){
            p1=which(GROUP %in% VP[i,1])
            p2=which(GROUP %in% VP[i,2])
            sd1=sd(THIS_PC[p1])
            mean1=mean(THIS_PC[p1])
            sd2=sd(THIS_PC[p2])
            mean2=mean(THIS_PC[p2])
            sd_lst1=c(sd_lst1,sd1)
            sd_lst2=c(sd_lst2,sd2)
            
            max_lst1=c(max_lst1,max(THIS_PC[p1]))
            min_lst1=c(min_lst1,min(THIS_PC[p1]))
            max_lst2=c(max_lst2,max(THIS_PC[p2]))
            min_lst2=c(min_lst2,min(THIS_PC[p2]))
            
            mean_lst1=c(mean_lst1,mean1)
            mean_lst2=c(mean_lst2,mean2)
            i=i+1}
    

        mean_com= apply(cbind(mean_lst1,mean_lst2),1,mean)    
   
      
        .x1_to_com=function(x1){
            if(x1<min(min_lst1)){x1=min(min_lst1)}
            if(x1>max(max_lst1)){x1=max(max_lst1)}
            x1=x1
            dlst1=c()
              
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst1[i]
                this_mean=mean_lst1[i]
                this_d=dnorm(x1,sd=this_sd,mean=this_mean)
                #######################   
                if(x1 > max_lst1[i] | x1 < min_lst1[i]){this_d=0}
                ##################
                if(is.na(this_d)){this_d=0}
                dlst1=c(dlst1,this_d)
                i=i+1} 
            
            if(sum(dlst1)==0){out=x1}else{out=sum(dlst1/sum(dlst1)*mean_com)}
            
            return(out)}
      
        .x2_to_com=function(x2){
            if(x2<min(min_lst2)){x2=min(min_lst2)}
            if(x2>max(max_lst2)){x2=max(max_lst2)}
            x2=x2
            dlst2=c()
            
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst2[i]
                this_mean=mean_lst2[i]
                this_d=dnorm(x2,sd=this_sd,mean=this_mean)
                #######################   
                if(x2 > max_lst2[i] | x2 < min_lst2[i]){this_d=0}
                ##################
                if(is.na(this_d)){this_d=0}
                dlst2=c(dlst2,this_d)
                i=i+1} 
           
            if(sum(dlst2)==0){out=x2}else{out=sum(dlst2/sum(dlst2)*mean_com)}
            
            return(out)}
         
         ########################
        
        lst1lst1=apply(as.matrix(DR[index1,THIS_DR]),1,.x1_to_com) 
        lst2lst2=apply(as.matrix(DR[index2,THIS_DR]),1,.x2_to_com)        
       
        OUT$adr[index1,THIS_DR]=lst1lst1
        OUT$adr[index2,THIS_DR]=lst2lst2
        #par(mfrow=c(1,2))
        #plot(DR[index1,THIS_DR],lst1lst1)
        #plot(DR[index2,THIS_DR],lst2lst2)
        
        lst1_mean=c()
        lst2_mean=c()
        i=1
        while(i<=nrow(VALID_PAIR)){
            this_pair=VALID_PAIR[i,]
            this_index1=which(GROUP %in% this_pair[1])
            this_index2=which(GROUP %in% this_pair[2])
            lst1_mean=c(lst1_mean,mean(OUT$adr[this_index1,THIS_DR]))
            lst2_mean=c(lst2_mean,mean(OUT$adr[this_index2,THIS_DR]))
            
            i=i+1}
        
        this_test=cor.test(lst1_mean,lst2_mean)#sum(dist_lst)
        this_ua_test=cor.test(mean_lst1, mean_lst2)
        
        this_cor=this_test$estimate
        this_pv=this_test$p.value
        
        this_un_cor=this_ua_test$estimate
        this_un_pv=this_ua_test$p.value
        
        ALL_COR=c(ALL_COR, this_cor)
        ALL_PV=c(ALL_PV, this_pv) 
        ALL_UCOR=c(ALL_UCOR,this_un_cor)   
        ALL_UPV=c(ALL_UPV,this_un_pv)
        print(THIS_DR)
        
        THIS_DR=THIS_DR+1}
    
    
    OUT$adjcor=ALL_COR
    OUT$adjpv=ALL_PV
    OUT$adjfdr=p.adjust(OUT$adjpv,method='fdr')
    OUT$cor=ALL_UCOR
    OUT$pv=ALL_UPV
    OUT$fdr=p.adjust(OUT$pv,method='fdr')
    print('Finished!!!')
    return(OUT)
   }




OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)
boxplot(OUT$cor,OUT$adjcor)


pbmc@dr$pca@cell.embeddings=DR
pbmc@dr$adjpca@cell.embeddings=OUT$adr


PCUSE=which(beerout$cor>0.9 & p.adjust(beerout$pv,method='fdr')<0.05) 

boxplot(OUT$cor[PCUSE],OUT$adjcor[PCUSE])

pbmc <- RunUMAP(object = pbmc, reduction.use='adjpca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
DimPlot(pbmc,reduction.use='umap',group.by='condition',pt.size=0.1)
DimPlot(pbmc,reduction.use='umap',group.by='map',pt.size=0.1)

pcpbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
DimPlot(pcpbmc,reduction.use='umap',group.by='condition',pt.size=0.1)
DimPlot(pcpbmc,reduction.use='umap',group.by='map',pt.size=0.1)


pbmc <- RunTSNE(object = pbmc, reduction.use='adjpca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
DimPlot(pbmc,reduction.use='tsne',group.by='condition',pt.size=0.1)
DimPlot(pbmc,reduction.use='tsne',group.by='map',pt.size=0.1)


#PCUSE=which(OUT$cor>0.9 & p.adjust(OUT$pv,method='fdr')<0.05) 
pcpbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
DimPlot(pcpbmc,reduction.use='tsne',group.by='condition',pt.size=0.1)
DimPlot(pcpbmc,reduction.use='tsne',group.by='map',pt.size=0.1)





#C1=OUT$cor
#C2=OUT$cor
#C3=OUT$cor
#C4=OUT$cor
#boxplot(C1,C2,C3,C4,OUT$ucor)


OUT$cor-OUT$ucor

pbmc=bastout$seurat



pbmc@dr$pca@cell.embeddings=OUT$adr


PCUSE=which(beerout$cor>0.8 & p.adjust(OUT$pv,method='fdr')<0.05) 
#PCUSE=which(OUT$ucor>0.8 & p.adjust(OUT$upv,method='fdr')<0.05)                          

                             #PCUSE=1:ncol(DR)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
pbmc <- RunUMAP(object = pbmc, reduction.use='oldpca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)

                             
DimPlot(pbmc,reduction.use='umap',group.by='condition',pt.size=0.1)
DimPlot(pbmc,reduction.use='umap',group.by='map',pt.size=0.1)

DimPlot(pbmc,reduction.use='umap',pt.size=0.1)



















load('TSNE.RData')
library('Seurat')
source('scRef.R')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]


USE=which(pbmc@meta.data$ori=='astrocytes_ependymal')

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

ref_tag=cbind(names(pbmc@ident), as.character(pbmc@meta.data$ori))    
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]
exp_sc_mat= exp_ref_mat[,USE]

getRanGene <- function(X){
    POS = which(X >0 )
    N=length(POS)/2
    KEEP = sample(x=POS, size=N )
    NEG = POS[which(!POS %in% KEEP)]
    X[NEG]=0
    return(X)
    }

set.seed(123)
sim_exp_sc_mat = apply(exp_sc_mat,2, getRanGene)


D1=sim_exp_sc_mat
D2=exp_ref_mat
colnames(D1)=paste0('sim_',colnames(D1))

#bastout=BAST(D1, D2, CNUM=10, PCNUM=50, FDR=1, COR=0, CPU=4, print_step=10)
bastout=BAST(D1, D2, CNUM=10, PCNUM=50, CPU=4, print_step=10)

DimPlot(bastout$seurat,reduction.use='umap',group.by='condition',pt.size=0.1)

plot(bastout$cor)
LABEL=c(rep('SIM_astrocytes_ependymal',ncol(D1)),as.character(ori_label[,2]))
bastout$seurat@meta.data$lab=LABEL
DimPlot(bastout$seurat,reduction.use='umap',group.by='lab',pt.size=0.1)

DimPlot(bastout$seurat,reduction.use='umap',group.by='condition',pt.size=0.1)





























D1X=.data2one(D1)
D2X=.data2one(D2)


G1=.getGroup(D1X,'D1',CNUM=10)
G2=.getGroup(D2X,'D2',CNUM=10)

VP_OUT=.getValidpair(D1, G1, D2, G2, CPU=4, method='kendall', print_step=10)
VP=VP_OUT$vp



colnames(D1)=paste0('sim_',colnames(D1))

EXP=.simple_combine(D1,D2)$combine
GROUP=c(G1,G2)
CONDITION=c(rep('D1',ncol(D1)),rep('D2',ncol(D2)))

library(Seurat)
pbmc=CreateSeuratObject(raw.data = EXP, min.cells = 0, min.genes = 0, project = "ALL")
pbmc@meta.data$group=GROUP
pbmc@meta.data$condition=CONDITION

MAP=rep('NA',length(GROUP))
MAP[which(GROUP %in% VP[,1])]='D1'
MAP[which(GROUP %in% VP[,2])]='D2'
pbmc@meta.data$map=MAP

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"))


PCNUM=50
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)


DR=pbmc@dr$pca@cell.embeddings
B1index=which(CONDITION=='D1')
B2index=which(CONDITION=='D2')


OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)


par(mfrow=c(1,2))
plot(OUT$cor,pch=16)
plot(-log(OUT$pv),pch=16)


pbmc@dr$oldpca=pbmc@dr$pca
pbmc@dr$pca@cell.embeddings=OUT$adr


PCUSE=which(p.adjust(OUT$pv,method='fdr')<0.05 & OUT$cor>0.8)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE)


LABEL=c(rep('SIM_astrocytes_ependymal',ncol(D1)),as.character(ori_label[,2]))
pbmc@meta.data$lab=LABEL

pdf('our_MAP.pdf',width=10,height=7)
DimPlot(object =pbmc, reduction.use = "umap", group.by = "map",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "umap", group.by = "condition",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "umap", group.by = "lab",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "umap",  pt.size = 0.1, do.return = TRUE)
dev.off()

##############
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scran", version = "3.8")
library(scran)

gene.counts1=D1
sce1 <- SingleCellExperiment(list(counts=gene.counts1))
sce1 <- normalize(sce1)

gene.counts2=D2
sce2 <- SingleCellExperiment(list(counts=gene.counts2))
sce2 <- normalize(sce2)

b1 <- sce1
b2 <- sce2
out <- fastMNN(b1, b2)
dim(out$corrected)

saveRDS(out,file='MNNout.RDS')

source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/BAST.R')

colnames(D1)=paste0('sim_',colnames(D1))

bastout=BAST(D1,D2,CNUM=10,PCNUM=50, FDR=0.05, COR=0.6)

saveRDS(bastout,file='BASTout.RDS')

TSNEPlot(bastout$seurat,group.by='condition')




pbmc=bastout$seurat
pbmc@dr$pca@cell.embeddings=out$corrected
rownames(pbmc@dr$pca@cell.embeddings)=rownames(pbmc@dr$oldpca@cell.embeddings)
colnames(pbmc@dr$pca@cell.embeddings)=colnames(pbmc@dr$oldpca@cell.embeddings)
PCUSE=1:50
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE)

TSNEPlot(pbmc,group.by='condition')
DimPlot(object =pbmc, reduction.use = "umap", group.by = "condition",  pt.size = 0.1, do.return = TRUE)
LABEL=c(rep('SIM_astrocytes_ependymal',ncol(D1)),as.character(ori_label[,2]))
pbmc@meta.data$lab=LABEL
DimPlot(object =pbmc, reduction.use = "umap", group.by = "condition",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "umap", group.by = "lab",  pt.size = 0.1, do.return = TRUE)





#PCUSE=1:50
#bastout$seurat <- RunUMAP(object = bastout$seurat, reduction.use='pca',dims.use = PCUSE)

bastout$seurat@meta.data$lab=LABEL
DimPlot(object =bastout$seurat, reduction.use = "umap", group.by = "lab",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =bastout$seurat, reduction.use = "umap", group.by = "condition",  pt.size = 0.1, do.return = TRUE)



