

EXP=pbmc.data
GENE=rownames(EXP)
ALL=c(as.character(LR[,1]),as.character(LR[,2]) )


expPermutate <- function(EXP, LRgene, TIME=10000,SEED=123){
  
    EXP=EXP
    SEED=SEED
    set.seed(SEED)
    ALL=LRgene
    GENE=rownames(EXP)  
    TIME=TIME
    permu_gene_index=which(GENE %in% ALL)

    MEAN=matrix(nrow=nrow(EXP[permu_gene_index,]),ncol=TIME)
    MEAN[which(is.na(MEAN))]=0
    rownames(MEAN)=rownames(EXP[permu_gene_index,])
    colnames(MEAN)=as.character(c(1:TIME))
    i=1
    while(i<=TIME){
        this_index=sample(c(1:ncol(EXP)),WINDOW)
        this_mean=apply(EXP[permu_gene_index,this_index],1,mean)
        MEAN[,i]=this_mean
        if(i%%100==1){print(i)}
        i=i+1
    }
    return(MEAN)}

#saveRDS(MEAN,file=paste0('MEAN.RDS' ))
