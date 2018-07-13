a=read.table('magic.csv',sep='\t',row.names=1,header=T,check.names=F)
a=t(a)
write.table(a,file='tmagic.csv',sep='\t',quote=F,row.names=T,col.names=T)
