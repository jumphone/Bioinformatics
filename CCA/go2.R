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



source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/scPA.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
CT=readRDS('CT.RDS')
MS=readRDS('MS.RDS')
CTX=readRDS('CTX.RDS')
MSX=readRDS('MSX.RDS')


CTG=.getGroup(CTX,'CT',CNUM=100)
MSG=.getGroup(MSX,'MS',CNUM=100)

VP_OUT=.getValidpair(CT, CTG, MS, MSG, CPU=4, method='kendall', do.plot=FALSE, print_step=10)
VP=VP_OUT$vp

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
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"))


PCNUM=50
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)

DR=pbmc@dr$pca@cell.embeddings
B1index=which(CONDITION=='CT')
B2index=which(CONDITION=='MS')



OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)

par(mfrow=c(1,2))
plot(OUT$cor,pch=16)
plot(-log(OUT$pv),pch=16)


pbmc@dr$oldpca=pbmc@dr$pca
pbmc@dr$pca@cell.embeddings=OUT$adr


PCUSE=which(p.adjust(OUT$pv,method='fdr')<0.05 & OUT$cor>0.8)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE)

pdf('MAP.pdf',width=10,height=7)
DimPlot(object =pbmc, reduction.use = "umap", group.by = "map",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "umap", group.by = "condition",  pt.size = 0.1, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "umap",  pt.size = 0.1, do.return = TRUE)
dev.off()








