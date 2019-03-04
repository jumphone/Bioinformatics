

library(Seurat)
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


.getGroup=function(X,TAG,CNUM=100){
    DR=X
    RANK=rank(DR,ties.method='random')
    CUTOFF=CNUM #round(max(RANK)/GNUM)
    GROUP=rep('NA',length(RANK))
    i=1
    j=1
    while(i<=length(RANK)){
        GROUP[which(RANK==i)]=paste0(TAG,'_',as.character(j))
        if(i%%CUTOFF==1){j=j+1;print(j)}
        i=i+1}
    return(GROUP)
}


