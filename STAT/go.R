a=read.table('data.txt',header=T,row.names=1,sep='\t')
b=a[which(!is.na(a[,1]*a[,2])),]
tmp=t.test(b[,2],b[,1],paired=T)
true_stat=tmp$statistic

c=c(b[,1],b[,2])
c=c[which(!is.na(c))]

END=100000
stat_list=c()

set.seed(12345)
i=1
while(i<=END){
  d1=sample(c,5)
  d2=sample(c,5) 
  this_stat=t.test(d2,d1,paired=T)$statistic
  stat_list=c(stat_list,this_stat)  
  if(i%%100==1){print(i)}
i=i+1}


PV=length(which(stat_list>true_stat))/END


#PV=0.04668

