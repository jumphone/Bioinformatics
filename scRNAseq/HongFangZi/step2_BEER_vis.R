source('BEER.R')
#####################
mybeer=readRDS(file='mybeer.RDS')





#pbmc <- mybeer$seurat
#PCUSE=mybeer$select
#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat

#umap=BEER.bbknn(pbmc, PCUSE, NB=5, NT=10)
#pbmc@reductions$umap@cell.embeddings=umap


pbmc=readRDS(file='pbmc.enh.RDS')


pdf('~/Downloads/HFZ3.pdf',width=10,height=10)
FeaturePlot(pbmc,features=c('FGF10','CGA','CYP19A1','GH2'),ncol=2)
FeaturePlot(pbmc,features=c('DIO2','LAIR2','CD105','CD44'),ncol=2)
FeaturePlot(pbmc,features=c('AIF1','CD14','CSF1R','CD53'),ncol=2)
FeaturePlot(pbmc,features=c('CD34','CDH5','ICAM1','PLVAP'),ncol=2)
FeaturePlot(pbmc,features=c('DKK1','IGFBP1','PBL','STAT3'),ncol=2)
FeaturePlot(pbmc,features=c('CD146','PDGFRB','CD31','NKH1'),ncol=2)
FeaturePlot(pbmc,features=c('PTPRC','CCL5','B2M'),ncol=2)
dev.off()
#saveRDS(pbmc, file='pbmc.enh.RDS')




