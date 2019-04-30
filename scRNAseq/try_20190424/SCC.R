
getBIN <- function(ONE, NUM=100){

    RANK=rank(ONE)
    LENGTH=length(ONE)
    WINDOW= round(LENGTH/(NUM))
    BIN=c()
    i=1
    while(i<=round(LENGTH/(WINDOW))){
        this_index=which((i-1)*WINDOW< RANK & i*WINDOW>=RANK)
        BIN=cbind(BIN,this_index)
        i=i+1
     }
    TAG=rep(NA,length(ONE))
    i=1
    while(i<=ncol(BIN)){
        TAG[BIN[,i]]=i
        i=i+1
        } 
    OUT=list(BIN=BIN,TAG=TAG,WINDOW=WINDOW)
    return(OUT)
    }


getMEAN <- function(EXP, LR, NUM=100,TIME=10000, SEED=123){
  
    LRgene=c(as.character(LR[,1]),as.character(LR[,2]))
    EXP=EXP
    SEED=SEED
    set.seed(SEED)
    ALL=LRgene
    GENE=rownames(EXP)  
    TIME=TIME
    permu_gene_index=which(GENE %in% ALL)
    LENGTH=ncol(EXP)
    WINDOW= round(LENGTH/(NUM))
  
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


getPMAT <- function(EXP, LR, BIN, MEAN ){

    LRgene=c(as.character(LR[,1]),as.character(LR[,2]))
    GENE=rownames(EXP)
    permu_gene_index=which(GENE %in% LRgene)
    EXP_LR=EXP[permu_gene_index,]
    TIME=ncol(MEAN)

    ECDF=c()
    j=1
    while(j<=nrow(MEAN)){
        ECDF=c(ECDF,ecdf(as.numeric(MEAN[j,])))
        j=j+1}

    PMAT = MEAN[,c(1:ncol(BIN))]*0
    i=1
    while(i<=ncol(BIN)){
        this_bin_index=BIN[,i]
  
        this_bin_mean_exp=apply(EXP_LR[,this_bin_index],1,mean)  
        this_p_list=c()
        j=1
        while(j<=length(this_bin_mean_exp)){
            this_p=-log(1+1/TIME-ECDF[[j]](this_bin_mean_exp[j]),10)
            this_p_list=c(this_p_list,this_p)
            j=j+1
            }
        PMAT[,i]=this_p_list
        print(i)
        i=i+1
        }
    return(PMAT)
    }



getCMAT <- function(EXP,LR,PMAT){
    
    GENE=rownames(EXP)
    CMAT=PMAT[c(1:ncol(PMAT)),]*0
    rownames(CMAT)=colnames(CMAT)
    rownames(CMAT)=paste0('L_',rownames(CMAT))
    colnames(CMAT)=paste0('R_',colnames(CMAT))

    i=1
    while(i<=nrow(LR)){

        this_l=as.character(LR[i,1])
        this_r=as.character(LR[i,2])
        if(this_l %in% GENE & this_r %in% GENE){
            this_l_index=which(rownames(PMAT)==this_l)
            this_r_index=which(rownames(PMAT)==this_r)
            this_l_bin_index=1
            while(this_l_bin_index<=nrow(CMAT)){
                this_r_bin_index=1
                while(this_r_bin_index<=ncol(CMAT)){
                    CMAT[this_l_bin_index,this_r_bin_index]=CMAT[this_l_bin_index,this_r_bin_index]+ 
                    PMAT[this_l_index,this_l_bin_index] - PMAT[this_r_index,this_l_bin_index] + PMAT[this_r_index,this_r_bin_index] - PMAT[this_l_index,this_r_bin_index]

                    this_r_bin_index=this_r_bin_index+1
                    }      
                this_l_bin_index=this_l_bin_index+1
                } 
             }
        if(i%%10==1){print(i)}
        i=i+1}

    CMAT=as.matrix(CMAT)
    return(CMAT)
    }


getPAIR <- function(CMAT){

    CUTOFF=0 
    TMP=CMAT
    TMP[which(TMP<CUTOFF)]=0
    PAIR=c()
    SCORE=c()
    i=1
    while(i<=ncol(TMP)){

        if(max(TMP[,i])>0){
        this_r=i
        j=1
        while(j<=nrow(TMP)){
            if(TMP[j,i]>0){
                this_l=j
                PAIR=cbind(PAIR,c(this_l,this_r))
                SCORE=c(SCORE,TMP[j,i])}
            j=j+1}
        }

        i=i+1}

    PAIR=t(PAIR)
    colnames(PAIR)=c('L','R')

    PAIR=PAIR[order(SCORE,decreasing=T),]
    SCORE=SCORE[order(SCORE,decreasing=T)]
    

    PAIR=PAIR[which(SCORE>0),]
    SCORE=SCORE[which(SCORE>0)]
    RANK=rank(-SCORE, ties.method = c( "min"))
    OUT=list(PAIR=PAIR,SCORE=SCORE,RANK=RANK)
    return(OUT)
    }


