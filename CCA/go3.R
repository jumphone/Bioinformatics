#Batch Alignment of Imbalanced Single-cell Data
#ISLET

source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/scPA.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


load('TSNE.RData')
library('Seurat')
source('scRef.R')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]


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
b2 <- sce1
out <- fastMNN(b1, b2)
dim(out$corrected)

saveRDS(out,file='MNNout.RDS')

source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/BAST.R')

bastout=BAST(D1,D2,CNUM=10,PCNUM=50)

saveRDS(bastout,file='BASTout.RDS')

