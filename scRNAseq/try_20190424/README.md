
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/SCC_LOGO.png" width="200">

### Single-Cell Communication (SCC) Toolkit 

Author: Feng Zhang
 
Date: 20190501
    
### Usage:

    library('Seurat')
    library(dplyr)
    library(Matrix)
    
    source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')
    
    pbmc=load('pbmc.RDS') # load Seurat Object
    
    #---- !!! Changed in Seurat 3.0 !!! ----
    pbmc.raw.data=getSeuratRAW(pbmc@raw.data, pbmc@scale.data)
    pbmc.data=as.matrix(pbmc@scale.data)
    used_gene=pbmc@var.genes  
    #---------------------------------------
    # For Seurat==3.0, please use:
    # pbmc.raw.data=getSeuratRAW(pbmc@assays$RNA@counts, pbmc@assays$RNA@scale.data)
    # pbmc.data=as.matrix(pbmc@assays$RNA@scale.data)
    # used_gene=VariableFeatures(object = pbmc)
    #---------------------------------------
        
    pbmc.raw.data=pbmc.raw.data[which(rownames(pbmc.raw.data) %in% used_gene),]
    pbmc.data=pbmc.data[which(rownames(pbmc.data) %in% used_gene),]
    
    #---- !!! Changed in Seurat 3.0 !!! ----
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')  
    #---------------------------------------
    # For Seurat==3.0, please use:
    # source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER_Seurat3.R')
    #---------------------------------------
    
    ONE=.data2one(pbmc.raw.data, used_gene, CPU=4, PCNUM=50, SEED=123,  PP=30)
    saveRDS(ONE,file='ONE.RDS')
    
    OUT=getBIN(ONE)
    BIN=OUT$BIN
    BINTAG=OUT$TAG
    saveRDS(BIN,file='BIN.RDS')
    saveRDS(BINTAG,file='BINTAG.RDS')
    
    pbmc@meta.data$bin=BINTAG
    pdf('1ID.pdf',width=12,height=10)
    DimPlot(pbmc,group.by='bin',reduction.use='umap',do.label=T)
    dev.off()
    
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/ID.png" width="200">

    LR=read.table('ReceptorLigand.txt.mouse',header=T,sep='\t')
    #https://github.com/jumphone/Bioinformatics/tree/master/scRNAseq/RecLig/
    
    EXP=pbmc.data
    
    MEAN=getMEAN(EXP, LR)
    saveRDS(MEAN,file='MEAN.RDS')
        
    PMAT=getPMAT(EXP, LR, BIN, MEAN)
    saveRDS(PMAT,file='PMAT.RDS')
    
    pdf('2GCOR.pdf',width=20,height=20)
    OUT=getPmatHEAT(PMAT,SHOW=T)
    dev.off()
    HEAT=OUT$HEAT
    DIST=OUT$DIST

    pdf('3CLUST.pdf',width=20,height=20)
    CLUST=getCLUST(HEAT, DIST, CCUT=0.7, SHOW=T)
    dev.off()
    HEAT=OUT$HEAT
    DIST=OUT$DIST

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/CLUST.png" width="200">

    MLR=getMLR(CLUST, LR, PMAT)
    LR=MLR[,c(1:2)]

    CMAT=getCMAT(EXP,LR,PMAT,BI=TRUE)
    saveRDS(CMAT,file='CMAT.RDS')
    
    pdf('4CMAT.pdf',width=15,height=13)
    library('gplots')
    heatmap.2(log(CMAT+1,10),scale=c("none"),dendrogram='both',Colv=T,Rowv=T,trace='none',
      col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
    heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',
      col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/CMAT.png" width="200">

    OUT=getPAIR(CMAT)
    PAIR=OUT$PAIR
    SCORE=OUT$SCORE
    RANK=OUT$RANK
    saveRDS(PAIR,file='PAIR.RDS')
   
    #---- !!! Changed in Seurat 3.0 !!! ----
    VEC=pbmc@dr$umap@cell.embeddings
    #--------------------------------------- 
    # For Seurat 3.0, please use:
    # VEC=pbmc@reductions$umap@cell.embeddings
    #---------------------------------------
    
    pdf('5CPlot_TOP200.pdf',width=12,height=10)
    CPlot(VEC,PAIR[1:200,],BINTAG)
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/CPlot.png" width="200">

    ORITAG=as.character(pbmc@ident)
    NET=getNET(PAIR, BINTAG,ORITAG )
    write.table(NET,file='NET.txt',sep='\t',row.names=F,col.names=T,quote=F)
       
    CN=getCN(NET)
    pdf('6DPlot.pdf',width=20,height=20)
    DP=DPlot(NET, CN, COL=3)
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/DPlot.png" width="200">

    SIG_INDEX=which(DP<0.05)
    SIG_PAIR=names(SIG_INDEX)
    
    pdf('7LPlot.pdf',width=20,height=20)
    RCN=trunc(sqrt(length(SIG_PAIR))+1)
    par(mfrow=c(RCN,RCN))
    i=1
    while(i<= length(SIG_PAIR) ){
        this_pair=SIG_PAIR[i]
        LT=unlist(strsplit(this_pair, "_to_"))[1]
        RT=unlist(strsplit(this_pair, "_to_"))[2]
        LP=LPlot(LT, RT, NET, PMAT,LR, MAIN=as.character(SIG_INDEX[i]),SEED=12345,PCUT=0.05)    
        colnames(LP)=paste0(c('Lexp','Rexp'),'_',c(LT,RT))
        write.table(LP,file=paste0(as.character(SIG_INDEX[i]),'.tsv'),row.names=T,col.names=T,sep='\t',quote=F)
        print(i)
        i=i+1}
    dev.off()
    
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/LPlot.png" width="200">


    TAG=groupTAG(BINTAG,LT="LGroup",RT='RGroup',LC=c(1,2,3),RC=c(4,5,6))
    TMP=getNET(PAIR, BINTAG, TAG )
    LPlot(LT='LGroup', RT='RGroup', TMP, PMAT,SEED=123)    
   
    

    
    
    
    
    
    
    
    
    
