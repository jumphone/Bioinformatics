a=read.table('test0.2.data.result.b.KM_100_SSS_TMP/tmp.6.k')
a=as.factor(t(a))
names(a)=names(EXP_cluster@ident)
EXP_cluster@ident=a
pdf('NEW.pdf',width=20,height=10)
TSNEPlot(object = EXP_cluster,do.label=T)
ALL=c(0:(length(table(EXP_cluster@ident))-1))
i=0
while(i<length(table(EXP_cluster@ident))){
print(i)
TSNEPlot(object = EXP_cluster,do.return = F,colors.use=c(rep('grey',length(table(EXP_cluster@ident))-1),'red'), plot.order=as.character(c(i,ALL[which(! ALL %in% i)])) )
i=i+1
}
dev.off()

