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
