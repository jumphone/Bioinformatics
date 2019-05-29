#Seurat 3.0

library(dplyr)
library(Seurat)



pbmc.data.1 <- Read10X(data.dir = "./CDC42_HET/")
pbmc.data.2 <- Read10X(data.dir = "./Small_Intestine/")

colnames(pbmc.data.1)=paste0('CDC42HET_',colnames(pbmc.data.1))
colnames(pbmc.data.2)=paste0('CDC42KO_',colnames(pbmc.data.2))

source('scRef.R')
DATA=.simple_combine(pbmc.data.1,pbmc.data.2)$combine

pbmc <- CreateSeuratObject(counts = DATA, project = "Intestine", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

pdf('QC1.pdf',width=15,height=5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 50)
pbmc@meta.data$batch=pbmc@meta.data$orig.ident



pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = c("percent.mt","nCount_RNA","batch"))
PCNUM=310
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=PCNUM)


#################
source('BEER_Seurat3.R')

EXP=as.matrix(pbmc@assays$RNA@counts)
D1=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='CDC42HET'),]) ) ]
D2=EXP[,which(colnames(EXP) %in%  rownames(pbmc@meta.data[which(pbmc@meta.data$batch=='CDC42KO'),]) ) ]


mybeer=BEER(D1, D2, CNUM=20, PCNUM=200, GN=5000, CPU=4, MTTAG="^mt-", REGBATCH=FALSE)
saveRDS(mybeer, file='mybeer.RDS')


plot(mybeer$cor, xlab='PCs', ylab="COR", pch=16)

npbmc <- mybeer$seurat
PCUSE <- which(mybeer$cor> quantile(mybeer$cor,0.2) & mybeer$fdr<0.05 )
npbmc <- RunUMAP(object = npbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)

PCUSE <- c(1:200)
allpbmc <- RunUMAP(object = npbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)


npbmc=readRDS('pbmc.RDS')

#npbmc@meta.data$group=as.character(npbmc@active.ident)
#npbmc@meta.data$group[which(npbmc@meta.data$group=='SmallIntestine')]='CDC42KO'

pdf('BEER.pdf',width=7,height=6)
DimPlot(allpbmc, reduction = "umap")
DimPlot(npbmc, reduction = "umap")
DimPlot(npbmc, reduction = "umap",group.by='map')
dev.off()

npbmc@meta.data$batch=npbmc@meta.data$orig.ident
saveRDS(npbmc,file='pbmc.RDS')

##############################
#source('scdemix.R')

#get_file(npbmc,TMP='./')
#agg_local(USED_CUTOFF=1, NUMBER=150, TMP='./')
#agg_cloud(npbmc,TMP="./")
#com_local(USED_CUTOFF=1,TMP='./')

#pbmc=readRDS('demix.RDS')


#pdf("GROUP.pdf",width=8,height=3)
#FeaturePlot(pbmc, cols=c("lightgrey",'red'), features = c('CDC42HET','SmallIntestine'))
#dev.off()

exp_ref1=read.table('MCA_INTEST.txt',header=T,row.names=1,sep='\t')
exp_ref2=read.table('GSE92332_intestin_mouse_ref.txt',header=T,row.names=1,sep='\t')



exp_sc=as.matrix(npbmc@assays$RNA@data)
source('scRef.R')
OUT1=SCREF(exp_sc, exp_ref1, CPU=4, min_cell=10,print_step=10)
#OUT1=OUT
npbmc@meta.data$mca=OUT1$tag2[,2]





pdf("MCA.pdf",width=10,height=5)
DimPlot(npbmc, group.by='mca',label=T)
dev.off()


pbmc=npbmc



pbmc@meta.data$type=pbmc@meta.data$mca
#pbmc@meta.data$type[which(pbmc@reductions$umap@cell.embeddings[,1] > -2.8)]='Epithelium'
pbmc@meta.data$type[which(pbmc@meta.data$mca %in% c('Columnar.epithelium_5',
"Epithelial_4","Epithelial.cell_17","Epithelium_9","Epithelium.of.small.intestinal.villi_13" ,                                          
"Epithelium.of.small.intestinal.villi_24","Epithelium.of.small.intestinal.villi_25",
"Epithelium.of.small.intestinal.villi_3","S.cell_16","S.cell_8","Stromal.cell_11","Macrophage_19" ,"Mast.cell_26"                                             
                               ))]='Epithelium'




pdf("MCA_EPI.pdf",width=10,height=5)
DimPlot(pbmc, group.by='type',label=T)
dev.off()

exp_sc=exp_sc[,which(pbmc@meta.data$type=='Epithelium')]
OUT2=SCREF(exp_sc, exp_ref2, CPU=4, min_cell=10,print_step=10)

pbmc@meta.data$all=pbmc@meta.data$type
pbmc@meta.data$all[which(pbmc@meta.data$type=='Epithelium')]=OUT2$tag2[,2]


pbmc@meta.data$all[which(pbmc@meta.data$all %in% c('Macrophage_19',
              'Macrophage_23','Macrophage_6'))]='Macrophage'

pbmc@meta.data$all[which(pbmc@meta.data$all %in% c('T.cell_10',
              'T.cell_12','T.cell_27','T.cell_7'))]='T.cell'


pbmc@meta.data$all[which(pbmc@meta.data$all %in% c('Dendrtic.cell_22'))]='Dendrtic.cell'
pbmc@meta.data$all[which(pbmc@meta.data$all %in% c('B.cell_2'))]='B.cell'

pdf("ALL.pdf",width=10,height=5)
DimPlot(pbmc, group.by='all',label=T)
dev.off()


CDC42HET=which(pbmc@active.ident=='CDC42HET')
#CDC42HET=rep(CDC42HET,pbmc@meta.data$CDC42HET[CDC42HET])


CDC42KO=which(pbmc@active.ident=='CDC42KO')
#CDC42KO=rep(CDC42KO,pbmc@meta.data$CDC42KO[CDC42KO])



write.table(sort(table(pbmc@meta.data$all[CDC42HET])),file='CDC42HET.txt',sep='\t',quote=F,col.names=F,row.names=F)
write.table(sort(table(pbmc@meta.data$all[CDC42KO])),file='SmallIntestine.txt',sep='\t',quote=F,col.names=F,row.names=F)

pdf("PIE.pdf",width=10,height=5)
pie(sort(table(pbmc@meta.data$all[CDC42HET])),main='DC42HET')
pie(sort(table(pbmc@meta.data$all[CDC42KO])),main='CDC42KO')
dev.off()
saveRDS(pbmc, file='ALL.RDS')
