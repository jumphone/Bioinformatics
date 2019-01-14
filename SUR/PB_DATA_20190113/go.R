SUR=read.table('SUR.txt',sep='\t',header=T)
used=which(!is.na(SUR[,2]))

library(survival)
library(survminer)

score=SUR[used,2]
this_sur=SUR[used,3]

plot(this_sur,score, xlab='OS (month)', ylab='Reponse Rate (%)', pch=16)
abline(h=0,col='red')
cor.test(this_sur, score, method='spearman')

highcut=0
lowcut=0

TYPE=rep('MED',length(used))
TYPE[which(score> highcut)]='High'
TYPE[which(score< lowcut)]='Low'
surtime=this_sur[which(TYPE!='MED')]
surevent=rep(1,length(which(TYPE!='MED')))
surtype=TYPE[which(TYPE!='MED')]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
ggsurvplot(fit, pval = TRUE)
surv_pvalue(fit)



####################################
