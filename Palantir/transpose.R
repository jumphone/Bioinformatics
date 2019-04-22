a=read.table('MGH54_mat.txt',header=T,row.names=1,sep='\t')
b=t(a)
write.table(b,file='t_MGH54_mat.txt',row.names=T,col.names=T,sep=',',quote=F)
