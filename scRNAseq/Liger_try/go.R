
setwd('C:/Users/cchmc/Desktop/BEER')
library(liger)
library(cowplot)



D1=readRDS('MGH36.RDS')
D2=readRDS('MGH53.RDS')
D3=readRDS('MGH54.RDS')
D4=readRDS('MGH60.RDS')
D5=readRDS('MGH93.RDS')
D6=readRDS('MGH97.RDS')

pbmc.data = list(D1=D1,D2=D2,D3=D3,D4=D4,D5=D5,D6=D6)

# Create liger object
a.pbmc <- createLiger(pbmc.data)
a.pbmc <- normalize(a.pbmc)
# Can pass different var.thresh values to each dataset if one seems to be contributing significantly
# more genes than the other
a.pbmc <- selectGenes(a.pbmc, var.thresh = c(0.3, 0.875), do.plot = F)
#s.var.genes <- #readRDS('~/Downloads/pbmc_alignment/var_genes.RDS')
#a.pbmc@var.genes <- s.var.genes
a.pbmc <- scaleNotCenter(a.pbmc)
#k.suggest <- suggestK(a.pbmc, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F,
#                     nrep = 5)
a.pbmc <- optimizeALS(a.pbmc, k=22, thresh = 5e-5, nrep = 3)

a.pbmc <- runTSNE(a.pbmc, use.raw = T)
p1 <- plotByDatasetAndCluster(a.pbmc, return.plots = T)
# Plot by dataset
#print(p1[[1]])

a.pbmc <- quantileAlignSNF(a.pbmc, resolution = 0.4, small.clust.thresh = 20)

a.pbmc <- runTSNE(a.pbmc)
p_a <- plotByDatasetAndCluster(a.pbmc, return.plots = T) 
# Modify plot output slightly
p_a[[1]] <- p_a[[1]] + theme_classic() + theme(legend.position = c(0.85, 0.15)) + 
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))

print(p_a[[1]])





