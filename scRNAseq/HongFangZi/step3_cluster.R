setwd('F:/HFZ')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

pbmc=readRDS('pbmc.enh.RDS')


tiff("./Marker/decidual.stomal.cell.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('DKK1','IGFBP1','PRL'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/epithelial.glandular.cell.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('STAT3','IGFBP1','STAT3'),ncol=2,pt.size=0.01)
dev.off()

tiff("./Marker/perivascular.cell.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('MCAM','PDGFRB','CD34','PECAM1'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/endothelial.cell.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD34','PLVAP','CDH5','ICAM1'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/NK.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('AMT','NCAM1','KLRK1'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/M.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD68','ITGAM','CD14','PTPRC'),ncol=2,pt.size=0.01)
dev.off()

tiff("./Marker/DC.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD14','CD52','CD83','CD86'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/T.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CCL5','CD3D','B2M','CXCR4'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/CTB.cytotrophoblast.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('KRT7','VGLL1','VGLL1'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/STB.syncytiotrophoblast.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CGA','CYP19A1','GH2'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/EVT.ExtravillousTrophoblast.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('HLA-G','DIO2','LAIR2'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/fibroblasts.MSC.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('ENG','VCAM1','SMN1','ITGA1'),ncol=2,pt.size=0.01)
dev.off()



tiff("./Marker/macrophages.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('AIF1','CD14','CD53','CSF1R'),ncol=2,pt.size=0.01)
dev.off()


tiff("./Marker/fetal.endothelial.tiff", width = 8, height= 9, units = 'in',res = 400)
FeaturePlot(pbmc, features=c('CD34','CDH5','ICAM1','PLVAP'),ncol=2,pt.size=0.01)
dev.off()



