#source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/TriTra.R')

#expmat=as.matrix(pbmc@assays$RNA@scale.data)
#UG=c('SOX11','SET')
#LG=c('VIM','SLC38A5')
#RG=c('NRL','TUBA1A')


TriTra <- function(expmat, UG,LG,RG,CUT=3){

    CUT=CUT
    expmat[which(expmat>CUT)]=CUT
    expmat[which(expmat< -CUT)]=-CUT

    UGE=expmat[which(rownames(expmat) %in% UG),]
    LGE=expmat[which(rownames(expmat) %in% LG),]
    RGE=expmat[which(rownames(expmat) %in% RG),]

    UGES=apply(UGE,2,mean)
    LGES=apply(LGE,2,mean)
    RGES=apply(RGE,2,mean)

    PU= UGES
    PL= LGES
    PR= RGES

    Y=PU-abs(PR-PL)
    X=PR-PL
    OUT=list()
    OUT$x=X
    OUT$y=Y
    return(OUT)
    }
