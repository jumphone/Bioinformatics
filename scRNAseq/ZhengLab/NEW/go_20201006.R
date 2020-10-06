

setwd('F:/Zhenglab/NewZhengZhang')

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
pbmc=readRDS(file='pbmc.final.RDS')



# Klf6

###############
library(ggplot2)
tiff(paste0("IMG/FX.Klf6.CDC42HET.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42HET'  )  ]
pbmc@meta.data$Klf6.EXP=pbmc@assays$RNA@data[which(rownames(pbmc)=='Klf6'),]
pbmc@meta.data$Klf6.EXP[which(!pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell'))]=0

p1=FeaturePlot(pbmc,features=c('Klf6.EXP'),cells=USED.CELL,order=TRUE,pt.size=2,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 6))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

dev.off()


tiff(paste0("IMG/FX.Klf6.CDC42KO.tiff"),width=5,height=4,units='in',res=600)
USED.CELL=colnames(pbmc)[which((!pbmc@meta.data$celltype %in% c('Immune.Cell')) & pbmc@meta.data$batch=='CDC42KO'  )  ]
pbmc@meta.data$Klf6.EXP=pbmc@assays$RNA@data[which(rownames(pbmc)=='Klf6'),]
pbmc@meta.data$Klf6.EXP[which(!pbmc@meta.data$celltype %in% c('TA.Cell','Stem.Cell'))]=0

p1=FeaturePlot(pbmc,features=c('Klf6.EXP'),cells=USED.CELL,order=TRUE,pt.size=2,combine = FALSE)
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'red'),  limits = c(0, 6))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

dev.off()
