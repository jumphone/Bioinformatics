library('Seurat')
pbmc=readRDS('GFPSN_merged.RDS')


#MCA
table(pbmc@meta.data$mca)
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='mca')


#Development
table(pbmc@meta.data$dev)
COLOR=heat.colors(n=length(table(pbmc@meta.data$dev))+2)
TSNEPlot(object = pbmc, colors.use=COLOR, group.by ='dev')

#Inj_bulk
table(pbmc@meta.data$injbulk)
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='injbulk')

#Inj_sc
table(pbmc@meta.data$injsc)
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='injsc')




