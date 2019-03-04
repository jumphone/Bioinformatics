##############################


library('pcaPP')
library(nloptr)
type_name=colnames(exp_ref_mat)

OUT=c()
i=1
while(i<=ncol(exp_sc_mat)){
########################
this_sc_exp=exp_sc_mat[,i]
   
eval_f <- function(x){
   pred_exp=exp_ref_mat %*% x
   real_exp=this_sc_exp
   y= - cor.fk(pred_exp, real_exp)
   return(y)
}
lb = rep(0.001,length(type_name))
ub = c()
for(tn in type_name){
   ub =c(ub , tmpR_max[which(names(tmpR_max)==tn)])
}
eval_g_eq <- function( x ) {
constr <- c( sum(x)-1 )
return(constr)
}
jub=jitter(ub)
x0=(jub)/sum(jub)
local_opts <- list( "algorithm" = "NLOPT_LN_PRAXIS","xtol_rel" = 1e-7 )
#local_opts <- list( "algorithm" = "NLOPT_LN_PRAXIS","xtol_rel" = 1.0e-7 )
opts <- list("algorithm"="NLOPT_LN_AUGLAG", "xtol_rel"= 1e-8,"maxeval" = 5000, "local_opts" = local_opts)
#opts <- list("algorithm"="NLOPT_GN_ISRES", "xtol_rel"=1.0e-8,"maxeval" = 10000)
res <- nloptr( x0=x0, eval_f=eval_f, lb=lb, ub=ub, eval_g_eq=eval_g_eq, opts=opts)
res

this_out=res$solution
########################
OUT=cbind(OUT,this_out)
print(i)
i=i+1}

OUT=t(OUT)
colnames(OUT)=type_name
rownames(OUT)=colnames(exp_sc_mat)
