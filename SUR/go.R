library(survival)
library(survminer)
score=
SUR=
TYPE=rep('MED',length(used))
TYPE[which(score> quantile(score,0.75))]='HighScore'
TYPE[which(score< quantile(score,0.25))]='LowScore'
surtime=SUR[which(TYPE!='MED')]
surevent=rep(1,length(which(TYPE!='MED')))
surtype=TYPE[which(TYPE!='MED')]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
ggsurvplot(fit, pval = TRUE)
surv_pvalue(fit)


library(survival)
library(survminer)
a=read.table('SUR.txt',header=T,sep='\t')
used=which(!is.na(a[,2]))
score=a[used,2]
SUR=a[used,4]
TYPE=rep('MED',length(used))
TYPE[which(score> 0 )]='HighScore'
TYPE[which(score< 0)]='LowScore'
surtime=SUR[which(TYPE!='MED')]
surevent=rep(1,length(which(TYPE!='MED')))
surtype=TYPE[which(TYPE!='MED')]
surtype=as.data.frame(surtype)
surv_object <- Surv(time = surtime, event = surevent)
fit <- survfit(surv_object ~ surtype, data=surtype)
ggsurvplot(fit, pval = TRUE)
surv_pvalue(fit)

