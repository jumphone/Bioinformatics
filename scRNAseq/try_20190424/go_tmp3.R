
R=3
N=100

RAND=runif(N)
r=runif(N)/5
centerX=rep(0,N)
centerY=rep(0,N)

#r = R * sqrt(RAND)
theta = RAND * 2 * pi

x = centerX + R * cos(theta)
y = centerY + R * sin(theta)

D1=cbind(x+r,y+r)


R=1
N=100

RAND=runif(N)
r=runif(N)/5
centerX=rep(0,N)
centerY=rep(0,N)

#r = R * sqrt(RAND)
theta = RAND * 2 * pi

x = centerX + R * cos(theta)
y = centerY + R * sin(theta)

D2=cbind(x+r,y+r)
#plot(x,y)


D=t(cbind(t(D1),t(D2)))
plot(D)

DD=dist(D)
DD=as.matrix(DD)

source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')

OUT=getPmatHEAT(DD,SHOW=T)
HEAT=OUT$HEAT
DIST=OUT$DIST
ORDER=HEAT$colInd

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')  
DIST=cor(DD,method='spearman')
OOO=.data2one(DD, colnames(DIST), CPU=4, PCNUM=50, SEED=123,  PP=30)
ORDER=order(OOO)



CLUST=getCLUST(ORDER, DIST, CCUT=0.7, SHOW=T)
TAG=CLUST[order(as.numeric(names(CLUST)))]
plot(D,col=TAG)

