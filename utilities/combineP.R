z2p =function(z){return(2*pnorm(abs(z), lower.tail = F))}

combineP=function(p){return(pchisq((sum(log(p))*-2), df=length(p)*2, lower.tail=F))}
