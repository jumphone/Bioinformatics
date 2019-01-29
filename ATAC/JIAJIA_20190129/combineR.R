
args <- commandArgs(TRUE)
CTRL <- args[1] #CTRL
KO <- args[2] #KO


a=read.table(CTRL,sep='\t',header=T,row.names=1)
b=read.table(KO,sep='\t',header=T,row.names=1)

source('scRef.R')


A=cbind(a[,4],a[,4])
B=cbind(b[,4],b[,4])
rownames(A)=rownames(a)
rownames(B)=rownames(b)
colnames(A)=c('Ctrl','Ctrl')
colnames(B)=c('KO','KO')

C=.simple_combine(A,B)

OUT=C$combine[,c(1,3)]
colnames(OUT)=c('Ctrl','KO')

write.table(OUT,file='combinedRPKM.txt',row.names=T,col.names=T,sep='\t',quote=F)
