library(Seurat)
ctrl.data <- read.table('run1663_5000_picard.bam.clean.bam.dge.txt',header=T,row.names=1, sep = "\t")
stim.data <- read.table('N709_6000_picard.bam.clean.bam.dge.txt' ,header=T,row.names=1,sep = "\t")

# Set up control object
ctrl <- CreateSeuratObject(raw.data = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl@meta.data$stim <- "CTRL"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)
# Set up stimulated object
stim <- CreateSeuratObject(raw.data = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim@meta.data$stim <- "STIM"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
stim <- NormalizeData(stim)
stim <- ScaleData(stim, display.progress = F)

# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))

immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:20)

immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:20)

pdf('TSNE.pdf',width=20,height=10)
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)
dev.off()

save.image('RData')
save(immune.combined, file = "Seurat.Robj")

#################
#################
#################


load('Seurat.Robj')

p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)







