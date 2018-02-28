source("https://bioconductor.org/biocLite.R")
biocLite("kstreet13/slingshot")

library(slingshot)
get_lineages(pcaX, clus, start.clus = 'm10')
