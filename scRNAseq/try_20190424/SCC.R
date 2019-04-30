
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
                    if(this_l_bin_index==this_r_bin_index){this_add=0}else{
                    this_add=PMAT[this_l_index,this_l_bin_index] - PMAT[this_r_index,this_l_bin_index] + PMAT[this_r_index,this_r_bin_index] - PMAT[this_l_index,this_r_bin_index]
                    }
                    CMAT[this_l_bin_index,this_r_bin_index]=CMAT[this_l_bin_index,this_r_bin_index]+ this_add
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
    
    CUTOFF=0
    PAIR=PAIR[which(SCORE>CUTOFF),]
    SCORE=SCORE[which(SCORE>CUTOFF)]
    RANK=rank(-SCORE, ties.method = c( "min"))
    OUT=list(PAIR=PAIR,SCORE=SCORE,RANK=RANK)
    return(OUT)
    }



CCPlot<-function(VEC, PAIR, BINTAG){
    
    library(cluster)
    BIN_FLAG=BINTAG
    plot(VEC,col='grey80',pch=16,cex=0.3,main=paste0('Cell Communication Plot (CCPlot)'))

    legend("topleft", legend=c("Ligand", "Recepter"),fill=c("green", "blue"))

    i=1
    while(i<=nrow(PAIR)){

        this_pair=PAIR[i,]
        this_l=which(BIN_FLAG==this_pair[1])
        this_r=which(BIN_FLAG==this_pair[2])
        this_l_vec=VEC[this_l,]
        this_r_vec=VEC[this_r,]

        
  
        start_point=pam(this_l_vec, 1)$medoids
        end_point= pam(this_r_vec, 1)$medoids
        size_ratio = (nrow(PAIR)-i+1)/nrow(PAIR)
        base_size=4

        transparent_ratio =150
        points(start_point[1],start_point[2],pch=16,cex=base_size*size_ratio, col=rgb(0, 255, 0, transparent_ratio, maxColorValue=255)  )
        points(end_point[1],end_point[2],pch=16,cex=base_size*size_ratio, col=rgb(0, 0, 255, transparent_ratio, maxColorValue=255) )

        points(this_l_vec,col='grey50',pch=16,cex=0.3)
        points(this_r_vec,col='grey50',pch=16,cex=0.3)

        text_col='red'
        text_cex=1
        text(x=start_point[1],y=start_point[2],label=as.character(this_pair[1]),pos=as.numeric(this_pair[1])%%4+1, col=text_col,cex=text_cex)  
        text(x=end_point[1],y=end_point[2],label=as.character(this_pair[2]),pos=as.numeric(this_pair[2])%%4+1, col=text_col,cex=text_cex)  
        segments(start_point[1], start_point[2], end_point[1],end_point[2],col='grey40',lty=3,lwd=1)
   
        i=i+1}

    }
