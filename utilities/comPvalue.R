source("https://bioconductor.org/biocLite.R")
biocLite("survcomp")
library('survcomp')
combine.test(c(0.1,0.1),method='fisher')
