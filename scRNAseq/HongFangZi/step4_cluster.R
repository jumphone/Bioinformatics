setwd('F:/HFZ')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

pbmc=readRDS('pbmc.enh.RDS')


VEC=pbmc@reductions$umap@cell.embeddings

# Here, we use K-means to do the clustering
N=150
set.seed(135)
K=kmeans(VEC,centers=N)

CLUST=K$cluster
pbmc@meta.data$clust=as.character(CLUST)

tiff("ID.tiff", width = 8, height= 8, units = 'in',res = 400)
DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)+NoLegend()
dev.off()

################################################################


pbmc@meta.data$decidua=rep(0,ncol(pbmc))
pbmc@meta.data$decidua[which(pbmc@meta.data$batch %in% c('decidua0117',
               'decidua0417.2','decidua508','decidua510','decidua514.2',                                      
                'decidua2018c','decidua2019c','decidua20190215','decidua20190420'         
                                                     ))]=1
tiff("decidua.tiff", width = 8, height= 8, units = 'in',res = 400)
FeaturePlot(pbmc, features='decidua',order=TRUE)
dev.off()


pbmc@meta.data$placenta=rep(1,ncol(pbmc))
pbmc@meta.data$placenta[which(pbmc@meta.data$batch %in% c('decidua0117',
               'decidua0417.2','decidua508','decidua510','decidua514.2',                                      
                'decidua2018c','decidua2019c','decidua20190215','decidua20190420'         
                                                     ))]=0
tiff("placenta.tiff", width = 8, height= 8, units = 'in',res = 400)
FeaturePlot(pbmc, features='placenta',order=TRUE)
dev.off()
#######################################



pbmc@meta.data$celltype=rep('epithelial.or.stomal',ncol(pbmc))


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(6,95,134,46,53,38,112,49,68,108,150))]='Endothelial'


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(88,144,113,79,97))]='NK'


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(126,88,145,45,135,102,80,100,90,77,59))]='T'


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(2,133,111))]='STB'


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(128,24,64,139,92,8,149,147,17))]='EVT'


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(63,141,86,85,61,51,36,16,42,44,75,125,103,115))]='CTB'


pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(117,9))]='perivascular'

pbmc@meta.data$celltype[which(pbmc@meta.data$clust %in% 
                              c(105,129,65,22,14))]=NA

DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.5,label=TRUE)


saveRDS(object = pbmc@meta.data,file='F:/HFZ/META.RDS')

####################



pbmc@meta.data$tag=rep('placenta',ncol(pbmc))
pbmc@meta.data$tag[which(pbmc@meta.data$batch %in% c('decidua0117',
               'decidua0417.2','decidua508','decidua510','decidua514.2',                                      
                'decidua2018c','decidua2019c','decidua20190215','decidua20190420'         
                                                     ))]='decidua'

TAB=table(pbmc@meta.data$celltype, pbmc@meta.data$tag)


TAB=table(pbmc@meta.data$celltype, pbmc@meta.data$batch)

.writeTable(TAB,PATH = 'F:/HFZ/TABLE.txt')



