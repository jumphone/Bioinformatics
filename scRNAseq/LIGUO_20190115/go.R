source('scRef.R')

T1=read.table('10_P19_cKO1.txt.uniq.txt',sep='\t',row.names=1,header=T)
T2=read.table('11_P19_cKO2.txt.uniq.txt',sep='\t',row.names=1,header=T)
T=.simple_combine(cbind(T1,T1),cbind(T2,T2))$combine
T=T[,c(1,3)]
N=read.table('Reference_expression.txt',sep='\t',row.names=1,header=TRUE)
REF=.simple_combine(T,N)$combine






