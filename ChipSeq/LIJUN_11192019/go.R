setwd('C:/Users/cchmc/Desktop/CUIQING')


a=read.csv('Latstumorcellspeaks.annotated.txt',header=T,sep='\t')

a.score=aggregate(a[,6],list(a[,16]),mean)



b=read.csv('mSCpeaks.annotated.txt',header=T,sep='\t')

b.score=aggregate(b[,6],list(b[,16]),mean)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')



A=cbind(a.score[,2],a.score[,2])
B=cbind(b.score[,2],b.score[,2])
rownames(A)=a.score[,1]
rownames(B)=b.score[,1]

DATA=.simple_combine(A,B)$combine[,c(1,3)]

colnames(DATA)=c('Lats','mSC')
.writeTable(DATA,PATH='combined.txt')




