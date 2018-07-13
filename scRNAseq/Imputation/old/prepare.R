exp_data=read.table('GSE70630_OG_processed_data_v2.txt.cleaned.txt',header=T,row.names=1,check.names=F)

exp_data=t(exp_data)
write.table(file='input.txt',exp_data,col.names=T,row.names=F,quote=F,sep='\t')
