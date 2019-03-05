# Batch Alignment of imbalanced Single-cell daTa (BAST)
# Author: Feng Zhang
# Date: Mar. 5, 2019


#library(Seurat)
#source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
#source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/BAST.R')

.data2one <- function(DATA, CPU=4, PCNUM=50){
    PCUSE=1:PCNUM
    print('Start')
    library(Seurat)
    print('Step1.Create Seurat Object...')
    DATA = CreateSeuratObject(raw.data = DATA, min.cells = 0, min.genes = 0, project = "DATA") 
    print('Step2.Normalize Data...')
    DATA <- NormalizeData(object = DATA, normalization.method = "LogNormalize", scale.factor = 10000)
    print('Step3.Scale Data...')
    DATA <- ScaleData(object = DATA, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
    print('Step4.PCA...')
    DATA <- RunPCA(object = DATA, pcs.compute=PCNUM, pc.genes = rownames(DATA@data), do.print = FALSE)
    print('Step5.tSNE...')
    DATA <- RunTSNE(object = DATA, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
    DR=DATA@dr$tsne@cell.embeddings
    print('Finished!!!')
    return(DR)
    }

.getGroup <- function(X,TAG,CNUM=100){
    DR=X
    RANK=rank(DR,ties.method='random')
    CUTOFF=CNUM 
    GROUP=rep('NA',length(RANK))
    i=1
    j=1
    while(i<=length(RANK)){
        GROUP[which(RANK==i)]=paste0(TAG,'_',as.character(j))
        if(i%%CUTOFF==1){j=j+1;print(j)}
        i=i+1}
    return(GROUP)
}


.getValidpair <- function(DATA1, GROUP1, DATA2, GROUP2, CPU=4, method='kendall', print_step=10){
    source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('Start')
    print('Step1.Generate Reference...')
    REF1=.generate_ref(DATA1, cbind(GROUP1, GROUP1), min_cell=1) 
    REF2=.generate_ref(DATA2, cbind(GROUP2, GROUP2), min_cell=1) 
    print('Step2.Calculate Correlation Coefficient...')
    out = .get_cor( REF1, REF2, method=method,CPU=CPU, print_step=print_step)
    print('Step3.Analyze Result...')
    tag1=.get_tag_max(out)
    tag2=.get_tag_max(t(out))
    V=c()
    i=1
    while(i<=nrow(tag1)){
        t1=tag1[i,1]
        t2=tag1[i,2]
        if(tag2[which(tag2[,1]==t2),2]==t1){V=c(V,i)}           
        i=i+1}
    VP=tag1[V,]
    C=c()
    t=1
    while(t<=nrow(VP)){
        this_c=out[which(rownames(out)==VP[t,2]),which(colnames(out)==VP[t,1])]
        C=c(C,this_c)
        t=t+1}
    #if(do.plot==TRUE){plot(C)}
    #VP=VP[which(C>=CUTOFF),]  
    print('Finished!!!')
    OUT=list()
    OUT$vp=VP
    OUT$cor=C
    return(OUT)
    }



.dr2adr <- function(DR, B1index, B2index, GROUP, VP, SEED=123){
    set.seed(SEED)
    library(dtw)
    library(MALDIquant)
    library(pcaPP)
    OUT=list()
    OUT$adr=DR
    VALID_PAIR=VP
    ALL_COR=c()   
    ALL_PV=c() 
    index1=B1index
    index2=B2index
  
    vindex1=which(GROUP %in% VP[,1])
    vindex2=which(GROUP %in% VP[,2])
    
    print('Start')
    THIS_DR=1
    while(THIS_DR<=ncol(DR)){
        THIS_PC = DR[,THIS_DR]
        M1=c()
        M2=c()
        maplst1=c()
        maplst2=c()
        i=1
        while(i<=nrow(VALID_PAIR)){
            this_pair=VALID_PAIR[i,]
            this_index1=which(GROUP %in% this_pair[1])
            this_index2=which(GROUP %in% this_pair[2])
            seq1=sort(THIS_PC[this_index1])
            seq2=sort(THIS_PC[this_index2])
            this_aln=dtw(seq1,seq2,keep=TRUE)
            maplst1=c(maplst1, seq1[this_aln$index1])
            maplst2=c(maplst2, seq2[this_aln$index2])
            
            i=i+1}
        comlst=cbind(maplst1,maplst2)
        compc=apply(comlst,1,mean)
        #compc=maplst2
        
        #plot(comlst[,1],compc)
        comlst1o=order(comlst[,1])
        .findlst1 <-function(x){y=sample(compc[which(comlst[,1]==x)],1);return(y)}
        vlst1lst1=  comlst[,1][comlst1o][match.closest(DR[index1,THIS_DR], comlst[,1][comlst1o])]
        lst1lst1=apply(as.matrix(vlst1lst1),1,.findlst1) 
                
        comlst2o=order(comlst[,2])
        .findlst2 <-function(x){y=sample(compc[which(comlst[,2]==x)],1);return(y)}
        vlst2lst2=  comlst[,2][comlst2o][match.closest(DR[index2,THIS_DR], comlst[,2][comlst2o])]
        lst2lst2=apply(as.matrix(vlst2lst2),1,.findlst2)
        #OO=order(THIS_PC[vindex2])
        #lst2lst2 = THIS_PC[vindex2][OO][match.closest(DR[index2,THIS_DR], THIS_PC[vindex2][OO])]
        #plot(lst2lst2,DR[index2,THIS_DR])
        
        OUT$adr[index1,THIS_DR]=lst1lst1
        OUT$adr[index2,THIS_DR]=lst2lst2
        
        lst1_mean=c()
        lst2_mean=c()
        i=1
        while(i<=nrow(VALID_PAIR)){
            this_pair=VALID_PAIR[i,]
            this_index1=which(GROUP %in% this_pair[1])
            this_index2=which(GROUP %in% this_pair[2])
            lst1_mean=c(lst1_mean,mean(OUT$adr[this_index1,THIS_DR]))
            lst2_mean=c(lst2_mean,mean(OUT$adr[this_index2,THIS_DR]))
            
            i=i+1}
        
        this_test=cor.test(lst1_mean,lst2_mean)#sum(dist_lst)
        this_cor=this_test$estimate
        this_pv=this_test$p.value
        
        ALL_COR=c(ALL_COR, this_cor)
        ALL_PV=c(ALL_PV, this_pv) 
        print(THIS_DR)
        THIS_DR=THIS_DR+1}
    
    #OUT$cor=ALL_COR
    OUT$cor=ALL_COR
    OUT$pv=ALL_PV
    print('Finished!!!')
    return(OUT)
    }



BAST <- function(D1, D2, CNUM=10, PCNUM=50, FDR=0.05, COR=0.8, CPU=4, print_step=10){
    RESULT=list()
    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    
    D1=D1
    D2=D2
    CNUM=CNUM
    PCNUM=PCNUM
    print_step=print_step
    
    print('Step1.Data to one-dimension...')
    D1X=.data2one(D1)
    D2X=.data2one(D2)

    G1=.getGroup(D1X,'D1',CNUM)
    G2=.getGroup(D2X,'D2',CNUM)
    
    print('Step2.Get Valid Pairs...')
    VP_OUT=.getValidpair(D1, G1, D2, G2, CPU, method='kendall', print_step)
    VP=VP_OUT$vp
    
    print('Step3.Combine Data...')
    EXP=.simple_combine(D1,D2)$combine
    GROUP=c(G1,G2)
    CONDITION=c(rep('D1',ncol(D1)),rep('D2',ncol(D2)))
    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='D1'
    MAP[which(GROUP %in% VP[,2])]='D2'
    
    print('Step4.Pre-process Combined Data...')
    pbmc=CreateSeuratObject(raw.data = EXP, min.cells = 0, min.genes = 0, project = "ALL")
    pbmc@meta.data$group=GROUP
    pbmc@meta.data$condition=CONDITION
    pbmc@meta.data$map=MAP
    pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)
    
    print('Step5.Subspace Alignment...')
    DR=pbmc@dr$pca@cell.embeddings 
    B1index=which(CONDITION=='D1')
    B2index=which(CONDITION=='D2')
    OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)
    pbmc@dr$oldpca=pbmc@dr$pca
    pbmc@dr$pca@cell.embeddings=OUT$adr
      
    print('Step6.UMAP & tSNE...')
    PCUSE=which(p.adjust(OUT$pv,method='fdr')<FDR & OUT$cor>COR)
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE)
    pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = PCUSE)
    ##########################
    
    RESULT$seurat=pbmc
    RESULT$vp=VP
    RESULT$d1x=D1X
    RESULT$d2x=D2X
    RESULT$g1=G1
    RESULT$g2=G2
    
    return(RESULT)
    }




