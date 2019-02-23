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

stat_list=stat_list[which(!is.na(stat_list))]

PV=length(which(stat_list>true_stat))/length(stat_list)

save.image(file='data.RData')
#PV=0.0466828

pdf('permutation_test.pdf',width=12,height=6)
#barplot(b[,1],b[,2])
par(mfrow=c(1,2))
M=c(mean(b[,1]),mean(b[,2]))
names(M)=c('primary','recurrence')
bp=barplot(M,ylim=c(0,20))

plot(density(stat_list),main='',xlab='Statistic',xlim=c(-4,4))
abline(v=true_stat,lwd=1.5,col='red')
text(true_stat,0.2,'p.value=0.046',pos=4,col='red')
#text(0,0.1,'100,000 permutations',pos=1,col='black')
dev.off()








