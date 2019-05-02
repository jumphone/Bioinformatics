###########################

# Author: Feng Zhang
# Date: 20190501

###########################

getSeuratRAW <- function(raw.data, scale.data){
    RAW=as.matrix(raw.data[,which(colnames(raw.data) %in% colnames(scale.data))])
    return(RAW)
    }

getBIN <- function(ONE, NUM=100){

    RANK=rank(ONE)
    LENGTH=length(ONE)
    WINDOW= trunc(LENGTH/(NUM))
    BIN=c()
    i=1
    while(i<=trunc(LENGTH/(WINDOW))){
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


BIN2BINTAG <- function(BIN, ONE){
    BINTAG=rep(NA,length(ONE))
    i=1
    while(i<=ncol(BIN)){
        BINTAG[BIN[,i]]=i
        i=i+1
        }
    return(BINTAG)
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
    WINDOW= trunc(LENGTH/(NUM))
  
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



getCMAT <- function(EXP, LR, PMAT, PRO=FALSE){
    
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
                    
                    l_bin_base = PMAT[this_l_index,this_l_bin_index] - PMAT[this_r_index,this_l_bin_index]
                    r_bin_base = PMAT[this_r_index,this_r_bin_index] - PMAT[this_l_index,this_r_bin_index]
                        
                    ######################  
                    if(PRO==FALSE){
                        this_add= l_bin_base + r_bin_base 
                    }else{
                        if(l_bin_base<=0 | r_bin_base<=0){
                            this_add=0}else{
                            this_add= l_bin_base + r_bin_base
                            }
                         }
                    ######################  
                        
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



CPlot<-function(VEC, PAIR, BINTAG){
    
    PAIR=matrix(PAIR, ncol=2)
    library(cluster)
    BIN_FLAG=BINTAG
    plot(VEC,col='grey80',pch=16,cex=0.3,main=paste0('Communication Plot (CPlot), N=',as.character(nrow(PAIR))))

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
        points(start_point[1],start_point[2],pch=16,cex=base_size*size_ratio+0.1, col=rgb(0, 255, 0, transparent_ratio, maxColorValue=255)  )
        points(end_point[1],end_point[2],pch=16,cex=base_size*size_ratio+0.1, col=rgb(0, 0, 255, transparent_ratio, maxColorValue=255) )

        points(this_l_vec,col='grey50',pch=16,cex=0.3)
        points(this_r_vec,col='grey50',pch=16,cex=0.3)

        text_col='red'
        text_cex=1
        text(x=start_point[1],y=start_point[2],label=as.character(this_pair[1]),pos=as.numeric(this_pair[1])%%4+1, col=text_col,cex=text_cex)  
        text(x=end_point[1],y=end_point[2],label=as.character(this_pair[2]),pos=as.numeric(this_pair[2])%%4+1, col=text_col,cex=text_cex)  
        segments(start_point[1], start_point[2], end_point[1],end_point[2],col='grey40',lty=3,lwd=1)
   
        i=i+1}

    }


getNET <- function(PAIR, BINTAG, ORITAG){
    
    source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    TT=table(ORITAG, BINTAG)
    tag=.get_tag_max(TT)

    PP=c()
    i=1
    while(i<=nrow(PAIR)){
        this_pair=PAIR[i,]
        this_tag1=tag[which(tag[,1]==as.character(this_pair[1])),2]
        this_tag2=tag[which(tag[,1]==as.character(this_pair[2])),2]
        PP=cbind(PP,c(this_tag1,this_tag2))
        i=i+1}
    PP=t(PP)
    
    
    
    SIZE=1-(1:nrow(PAIR))/(nrow(PAIR)+1)
    
    NET=cbind(PAIR,PP,SIZE)
    NET=as.matrix(NET)
    colnames(NET)=c('L','R','LT','RT','RANK')
    return(NET)
    }

getCN <- function(NET){
    OUTPUT=NET
    CN=paste0(OUTPUT[,3],'_to_',OUTPUT[,4])
    #par(mar=c(5,15,5,15))
    return(sort(table(CN),decreasing=TRUE))
    }

DPlot <- function(NET, CN, CUTOFF=3, PCUT=0.05, COL=2,PLOT=TRUE){   
    CCLR=names(CN[which(CN>=CUTOFF)])
    if(PLOT==TRUE){
        par(mfrow=c(trunc((length(CCLR)+1)/COL),COL))}
    TOT=paste0(as.character(NET[,3]),'_to_',as.character(NET[,4]))
    ALLP=c()
    i=1
    while(i<=length(CCLR)){
        this_cclr=CCLR[i]
        this_p=ks.test(which(TOT==this_cclr),jitter(1:length(TOT)),alternative='greater')$p.value
        if(this_p<PCUT){CM='red'}else{CM='black'}
        ALLP=c(ALLP,this_p)
        this_p=signif(this_p, digits = 2)
        this_p=format(this_p, scientific = T)
        if(PLOT==TRUE){            
            plot(col.main=CM, main=paste0(as.character(i),': ',this_cclr,'; KS p-value=',this_p),x=which(TOT==this_cclr),y=rep(1,length(which(TOT==this_cclr))),type='h',ylim=c(0,1),xlim=c(0,length(TOT)),xlab='RANK',ylab='',col='green3') 
        }
        i=i+1}
    names(ALLP)=CCLR
    return(ALLP)
    }

LPlot <- function(LT,RT,NET,PMAT,MAIN='',SEED=123){
    

    set.seed(SEED)
    OUTPUT=NET
    GENE=rownames(PMAT )
    VP=OUTPUT[which(OUTPUT[,3]==LT & OUTPUT[,4]==RT),]
    if(length(VP)>5){
    this_l_exp=apply(PMAT[,as.numeric(VP[,1])],1,mean)
    this_r_exp=apply(PMAT[,as.numeric(VP[,2])],1,mean)
    }else{this_l_exp=PMAT[,as.numeric(VP[1])];this_r_exp=PMAT[,as.numeric(VP[2])] }
    
    tag_list=c()
    #out_list=c()
    l_list=c()
    r_list=c()

    i=1
    while(i<=nrow(LR)){

        this_l=LR[i,1]
        this_r=LR[i,2]
        this_tag=paste0(this_l,"_",this_r)
        if(this_l %in% GENE & this_r %in% GENE){
            #this_out=this_l_exp[which(names(this_l_exp)==this_l)]+this_r_exp[which(names(this_r_exp)==this_r)]
            this_l_out=this_l_exp[which(names(this_l_exp)==this_l)]
            this_r_out=this_r_exp[which(names(this_r_exp)==this_r)]
            tag_list=c(tag_list,this_tag)
            #out_list=c(out_list,this_out)
            l_list=c(l_list,this_l_out)
            r_list=c(r_list,this_r_out)
            }  
   
        i=i+1
        }
    names(l_list)=paste0(tag_list,'_L')
    names(r_list)=paste0(tag_list,'_R')
    
    XLIM=c(0,max(as.numeric(r_list))+0.2)
    YLIM=c(0,max(as.numeric(l_list))+0.2)
    
    
    plot(main=paste0(MAIN,' Expression range: 0~',max(round(as.numeric(PMAT)))), r_list, l_list, pch=16, xlab=paste0('EXP of Receptor in ',RT),ylab=paste0('EXP of Ligend in ',LT) ,xlim=XLIM,ylim=YLIM)
    text(r_list, l_list, label=tag_list,pos=sample(c(1,2,3,4),length(l_list),replace = TRUE))
    OUT=cbind(l_list,r_list)
    rownames(OUT)=tag_list
    colnames(OUT)=c('Lexp','Rexp')
    return(OUT)
}




groupTAG <- function(BINTAG,LT,RT,LC,RC){
    LT=LT
    RT=RT
    ORITAG=rep('NA',length(BINTAG))
    ORITAG[which(BINTAG %in% LC)]=LT
    ORITAG[which(BINTAG %in% RC)]=RT
    return(ORITAG)  
}





