# Single-Cell Communication (SCC) Toolkit 

## Author: Feng Zhang
 
## Date: 20190501
    
        
    library('Seurat')
    library(dplyr)
    library(Matrix)
    
    source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/scRNAseq/try_20190424/SCC.R')
    
    pbmc=load('pbmc.RDS') # load Seurat Object
    
    pbmc.raw.data=getSeuratRAW(pbmc@raw.data,pbmc@scale.data)
    pbmc.data=as.matrix(pbmc@scale.data)
    used_gene=pbmc@var.genes   
    
    # For Seurat==3.0, please use:
    # pbmc.raw.data=getSeuratRAW(pbmc@assays$RNA@counts, pbmc@assays$RNA@scale.data)
    # pbmc.data=as.matrix(pbmc@assays$RNA@scale.data)
    # used_gene=VariableFeatures(object = pbmc)
    
    pbmc.raw.data=pbmc.raw.data[which(rownames(pbmc.raw.data) %in% used_gene),]
    pbmc.data=pbmc.data[which(rownames(pbmc.data) %in% used_gene),]
    
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')  
    
    # For Seurat==3.0, please use:
    # source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER_Seurat3.R')
    
    ONE=.data2one(pbmc.raw.data, used_gene, CPU=4, PCNUM=50, SEED=123,  PP=30)
    saveRDS(ONE,file='ONE.RDS')
    
    OUT=getBIN(ONE)
    BIN=OUT$BIN
    BINTAG=OUT$TAG
    saveRDS(BIN,file='BIN.RDS')
    saveRDS(BINTAG,file='BINTAG.RDS')
    
    pbmc@meta.data$bin=BINTAG
    pdf('1ID.pdf',width=12,height=10)
    DimPlot(pbmc,group.by='bin',reduction.use='tsne',do.label=T)
    dev.off()
    
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/ID.png" width="200">

    LR=read.table('RL_mouse.txt',header=T,sep='\t')
    EXP=pbmc.data
    
    MEAN=getMEAN(EXP, LR)
    saveRDS(MEAN,file='MEAN.RDS')
        
    PMAT=getPMAT(EXP, LR, BIN, MEAN)
    saveRDS(PMAT,file='PMAT.RDS')
    
    CMAT=getCMAT(EXP,LR,PMAT)
    saveRDS(CMAT,file='CMAT.RDS')
    
    pdf('2HEAT.pdf',width=15,height=13)
    library('gplots')
    heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',
      col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/HEAT.png" width="200">

    OUT=getPAIR(CMAT)
    PAIR=OUT$PAIR
    SCORE=OUT$SCORE
    RANK=OUT$RANK
    saveRDS(PAIR,file='PAIR.RDS')
   
    VEC=pbmc@dr$tsne@cell.embeddings
    
    # For Seurat 3.0, please use:
    # VEC=pbmc@reductions$tsne@cell.embeddings 
    
    pdf('3CPlot.pdf',width=12,height=10)
    CPlot(VEC,PAIR[1:200,],BINTAG)
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/CPlot.png" width="200">

    ORITAG=as.character(pbmc@ident)
    NET=getNET(PAIR, BINTAG,ORITAG )
    write.table(NET,file='NET.txt',sep='\t',row.names=F,col.names=T,quote=F)
       
    CN=getCN(NET)
    pdf('4DPlot.pdf',width=20,height=20)
    DP=DPlot(NET, CN, COL=3)
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/DPlot.png" width="200">

    SIG_INDEX=which(DP<0.05)
    SIG_PAIR=names(SIG_INDEX)
    
    pdf('5LPlot.pdf',width=20,height=20)
    RCN=trunc(sqrt(length(SIG_PAIR))+1)
    par(mfrow=c(RCN,RCN))
    i=1
    while(i<= length(SIG_PAIR) ){
        this_pair=SIG_PAIR[i]
        LT=unlist(strsplit(this_pair, "_to_"))[1]
        RT=unlist(strsplit(this_pair, "_to_"))[2]
        LP=LPlot(LT, RT, NET, PMAT,SEED=123)    
        colnames(LP)=paste0(c('Lexp','Rexp'),'_',c(LT,RT))
        write.table(LP,file=paste0(as.character(SIG_INDEX[i]),'.tsv'),row.names=T,col.names=T,sep='\t',quote=F)
        print(i)
        i=i+1}
    dev.off()
    
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/LPlot.png" width="200">
    
    
    
    
    
    
