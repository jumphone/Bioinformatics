
source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/scPA.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
CT=readRDS('CT.RDS')
MS=readRDS('MS.RDS')
CTX=readRDS('CTX.RDS')
MSX=readRDS('MSX.RDS')

CTG=.getGroup(CTX,'CT',CNUM=100)
MSG=.getGroup(MSX,'MS',CNUM=100)

VP=.getValidpair(CT, CTG, MS, MSG, CPU=4, method='kendall', do.plot=FALSE, print_step=10)




