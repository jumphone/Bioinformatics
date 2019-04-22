a=read.table('MGH54_mat.txt',header=T,row.names=1,sep='\t')
b=t(a)

varC=apply(b,2,var)


write.table(b[,which(varC>0.5)],file='t_MGH54_mat.txt',row.names=T,col.names=T,sep=',',quote=F)
