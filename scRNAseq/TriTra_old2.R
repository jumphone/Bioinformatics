#source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/TriTra.R')

#expmat=as.matrix(pbmc@assays$RNA@scale.data)
#UG=c('SOX11','SET')
#LG=c('VIM','SLC38A5')
#RG=c('NRL','TUBA1A')

#TriTra <- function(expmat, UG,LG,RG, delta=0.2){

 #   CUT=CUT
#    expmat[which(expmat>CUT)]=CUT
#    expmat[which(expmat< -CUT)]=-CUT

#    UGE=expmat[which(rownames(expmat) %in% UG),]
#    LGE=expmat[which(rownames(expmat) %in% LG),]
#    RGE=expmat[which(rownames(expmat) %in% RG),]

#    UGES=apply(UGE,2,mean)
#    LGES=apply(LGE,2,mean)
#    RGES=apply(RGE,2,mean)

#    PU= UGES
#    PL= LGES
#    PR= RGES

#    Y=PU-abs(PR-PL)
#    X=PR-PL
#    OUT=list()
#    OUT$x=X
#   OUT$y=Y
#    return(OUT)
#    }



TriTra <- function(expmat, UG,LG,RG,delta=0.2){
    ##################
    library(plot3D)
    ##################
    
    expmat=as.matrix(expmat)
    delta=delta
    
    #CUT=CUT
    #expmat[which(expmat>CUT)]=CUT
    #expmat[which(expmat< -CUT)]=-CUT

    UGE=expmat[which(rownames(expmat) %in% UG),]
    LGE=expmat[which(rownames(expmat) %in% LG),]
    RGE=expmat[which(rownames(expmat) %in% RG),]

    UGES=apply(UGE,2,mean)
    LGES=apply(LGE,2,mean)
    RGES=apply(RGE,2,mean)
  
    set.seed(123)
    RUGES=rank(UGES,ties.method = 'min')
    RLGES=rank(LGES,ties.method = 'min')
    RRGES=rank(RGES,ties.method = 'min')
    RUGES[which(UGES==0)]=0
    RLGES[which(LGES==0)]=0
    RRGES[which(RGES==0)]=0
    RUGES=RUGES/max(RUGES)
    RLGES=RLGES/max(RLGES)
    RRGES=RRGES/max(RRGES)
  
  
    SMAT=cbind(RUGES,RLGES,RRGES)
     
    SUGES=RUGES-apply(SMAT[,c(2,3)],1,max)
    SUGES[which(SUGES<0)]=0
    
    SLGES=RLGES-apply(SMAT[,c(1,3)],1,max)
    SLGES[which(SLGES<0)]=0
  
    SRGES=RRGES-apply(SMAT[,c(1,2)],1,max)
    SRGES[which(SRGES<0)]=0
  
    #DEL=quantile(c(SUGES,SLGES,SRGES),0.9)
    
    PU= RUGES * (delta+SUGES)
    PL= RLGES * (delta+SLGES)
    PR= RRGES * (delta+SRGES)

    
    pdf('./TT.tmp.pdf',width=5,height=5)
    pmat <- scatter3D(x=PL, y=PR,z=PU,  cex = 1, colkey = FALSE,theta = 130, phi = 30)
    dev.off()
    D2=trans3d(x=PL, y=PR,z=PU,pmat)
    
  
    OUT=list()
    OUT$x=D2$x
    OUT$y=D2$y
    return(OUT)
    }
