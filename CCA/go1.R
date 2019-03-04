source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/scPA.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

CT=readRDS('CT.RDS')
MS=readRDS('MS.RDS')

CTX=.data2one(CT)
MSX=.data2one(MS)


saveRDS(CT,file='CT.RDS')
saveRDS(MS,file='MS.RDS')
saveRDS(CTX,file='CTX.RDS')
saveRDS(MSX,file='MSX.RDS')


CTG=.getGroup(CTX,'CT',CNUM=100)
MSG=.getGroup(MSX,'MS',CNUM=100)

VP=.getValidpair(CT, CTG, MS, MSG, CPU=4, method='kendall', do.plot=FALSE, print_step=10)
  

EXP=.simple_combine(CT,MS)$combine
GROUP=c(CTG,MSG)
CONDITION=c(rep('CT',ncol(CT)),rep('MS',ncol(MS)))

library(Seurat)
pbmc=CreateSeuratObject(raw.data = EXP, min.cells = 0, min.genes = 0, project = "ALL")
pbmc@meta.data$group=GROUP
pbmc@meta.data$condition=CONDITION

MAP=rep('NA',length(GROUP))
MAP[which(GROUP %in% VP[,1])]='CT'
MAP[which(GROUP %in% VP[,2])]='MS'
pbmc@meta.data$map=MAP


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"))

PCNUM=100
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)


DR=pbmc@dr$pca@cell.embeddings
B1index=which(CONDITION=='CT')
B2index=which(CONDITION=='MS')


OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)



pbmc@dr$alnpca=pbmc@dr$pca
pbmc@dr$alnpca@key='APC'
pbmc@dr$alnpca@cell.embeddings=OUT$adr

#RATIO=0.8
#PCUSE=which(p.adjust(OUT$pv[2,],method='fdr')<0.05 & OUT$coef[2,] >0)

PCUSE=which(OUT$coef[2,] <1/RATIO & OUT$coef[2,]>RATIO & OUT$pv[2,] <0.05)

pbmc <- RunTSNE(object = pbmc, reduction.use='alnpca',dims.use = PCUSE, do.fast = TRUE)

DimPlot(object =pbmc, reduction.use = "tsne", group.by = "map",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "tsne", group.by = "condition",  pt.size = 0.1, do.return = TRUE)




