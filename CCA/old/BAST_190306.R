# Batch Alignment of imbalanced Single-cell daTa (BAST)
# Author: Feng Zhang
# Date: Mar. 5, 2019


#library(Seurat)
#library(dtw)
#library(MALDIquant)
#library(pcaPP)
#source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
#source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/BAST.R')

.data2one <- function(DATA, GENES, CPU=4, PCNUM=50){
    PCUSE=1:PCNUM
    print('Start')
    library(Seurat)
    print('Step1.Create Seurat Object...')
    DATA = CreateSeuratObject(raw.data = DATA, min.cells = 0, min.genes = 0, project = "DATA") 
    print('Step2.Normalize Data...')
    DATA <- NormalizeData(object = DATA, normalization.method = "LogNormalize", scale.factor = 10000)
    print('Step3.Scale Data...')
    DATA <- ScaleData(object = DATA, genes.use =GENES, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
    print('Step4.PCA...')
    DATA <- RunPCA(object = DATA, pcs.compute=PCNUM, pc.genes =GENES, do.print = FALSE)
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
        #if(i%%CUTOFF==1){j=j+1;print(j)}
        if(i%%CUTOFF==1){j=j+1}
        i=i+1}
    print(j-1)
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
       
        ########################
        sd_lst1=c()
        mean_lst1=c()
        max_lst1=max(THIS_PC[vindex1])
        min_lst1=min(THIS_PC[vindex1])
        sd_lst2=c()
        mean_lst2=c()  
        max_lst2=max(THIS_PC[vindex2])
        min_lst2=min(THIS_PC[vindex2])
        i=1
        while(i<=nrow(VP)){
            p1=which(GROUP %in% VP[i,1])
            p2=which(GROUP %in% VP[i,2])
            sd1=sd(THIS_PC[p1])
            mean1=mean(THIS_PC[p1])
            sd2=sd(THIS_PC[p2])
            mean2=mean(THIS_PC[p2])
            sd_lst1=c(sd_lst1,sd1)
            sd_lst2=c(sd_lst2,sd2)
            mean_lst1=c(mean_lst1,mean1)
            mean_lst2=c(mean_lst2,mean2)
            i=i+1}
    

        mean_com= apply(cbind(mean_lst1,mean_lst2),1,mean)    
   
      
        .x1_to_com=function(x1){
            #if(x1 <=min_lst1){x1=min_lst1}
            #if(x1 >=max_lst1){x1=max_lst1}
                
            x1=x1
            dlst1=c()
            value1=c()
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst1[i]
                this_mean=mean_lst1[i]
                this_d=dnorm(x1,sd=this_sd,mean=this_mean)
                #######################   
                this_v=(x1-this_mean)+mean_com[i]
                value1=c(value1,this_v)
                ##################
                if(is.na(this_d)){this_d=0}
                dlst1=c(dlst1,this_d)
                i=i+1} 
            
            out=sum(dlst1/sum(dlst1)*value1)
            return(out)}
      
        .x2_to_com=function(x2){
            #if(x2 <=min_lst2){x2=min_lst2}
            #if(x2 >=max_lst2){x2=max_lst2}
                
            x2=x2
            dlst2=c()
            value2=c()
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst2[i]
                this_mean=mean_lst2[i]
                this_d=dnorm(x2,sd=this_sd,mean=this_mean)
                #######################   
                this_v=(x2-this_mean)+mean_com[i]
                value2=c(value2,this_v)
                ##################
                if(is.na(this_d)){this_d=0}
                dlst2=c(dlst2,this_d)
                i=i+1} 
            out=sum(dlst2/sum(dlst2)*value2)
            return(out)}
         
         ########################
         #tmp=c(-10:20)
         #tmp1=apply(as.matrix(tmp),1,.x1_to_com)
         #tmp2=apply(as.matrix(tmp),1,.x2_to_com) 
         #plot(tmp1,tmp2)
        
        lst1lst1=apply(as.matrix(DR[index1,THIS_DR]),1,.x1_to_com) 
        lst2lst2=apply(as.matrix(DR[index2,THIS_DR]),1,.x2_to_com)        
        
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
        
        #this_fit=lm(lst2_mean~lst1_mean)
        #this_coef=this_fit$coefficients
        #OUT$adr[index1,THIS_DR]=this_coef[1]+this_coef[2]*OUT$adr[index1,THIS_DR]
        
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



BAST <- function(D1, D2, CNUM=10, PCNUM=50, CPU=4, print_step=10){
    RESULT=list()
    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BAST Start ')
    D1=D1
    D2=D2
    CNUM=CNUM
    PCNUM=PCNUM
    print_step=print_step
      
    print('######################################')
    print('MainStep1.Combine Data...')
    print('######################################')
    EXP=.simple_combine(D1,D2)$combine
   
    pbmc=CreateSeuratObject(raw.data = EXP, min.cells = 0, min.genes = 0, project = "ALL")
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    #length(x = pbmc@var.genes)
    pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
    pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)
    
    print('######################################')
    print('MainStep2.Data to one-dimension...')
    print('######################################')
    D1X=.data2one(D1, pbmc@var.genes, CPU, PCNUM)
    D2X=.data2one(D2, pbmc@var.genes, CPU, PCNUM)
    G1=.getGroup(D1X,'D1',CNUM)
    G2=.getGroup(D2X,'D2',CNUM)
    GROUP=c(G1,G2)
    CONDITION=c(rep('D1',ncol(D1)),rep('D2',ncol(D2)))
    pbmc@meta.data$group=GROUP
    pbmc@meta.data$condition=CONDITION
    
    
    print('######################################')
    print('MainStep3.Get Valid Pairs...')
    print('######################################')
    VP_OUT=.getValidpair(D1, G1, D2, G2, CPU, method='kendall', print_step)
    #VP_OUT=.getValidpair(D1, G1, D2, G2, 4, 'kendall', 10)
    VP=VP_OUT$vp
    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='D1'
    MAP[which(GROUP %in% VP[,2])]='D2'
    pbmc@meta.data$map=MAP
    
    print('######################################')
    print('MainStep4.Subspace Alignment...')
    print('######################################')
    DR=pbmc@dr$pca@cell.embeddings 
    B1index=which(CONDITION=='D1')
    B2index=which(CONDITION=='D2')
    OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)
    pbmc@dr$oldpca=pbmc@dr$pca
    pbmc@dr$pca@cell.embeddings=OUT$adr
    
    #print('######################################')
    #print('MainStep5.UMAP & tSNE...')
    #print('######################################')
    #PCUSE=1:PCNUM #which(p.adjust(OUT$pv,method='fdr')<FDR & OUT$cor>COR )
    #pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE)
    #pbmc <- RunUMAP(object = pbmc, reduction.use='oldpca',dims.use = PCUSE)
    #pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = PCUSE)
    #DimPlot(pbmc,reduction.use='umap',group.by='condition',pt.size=0.1)
    #DimPlot(pbmc,reduction.use='umap',group.by='map',pt.size=0.1)
    #DimPlot(pbmc,reduction.use='umap',pt.size=0.1)
    ########################## 
    RESULT$seurat=pbmc
    RESULT$vp=VP
    RESULT$d1x=D1X
    RESULT$d2x=D2X
    RESULT$g1=G1
    RESULT$g2=G2
    RESULT$cor=OUT$cor
    RESULT$pv=OUT$pv
    RESULT$fdr=p.adjust(OUT$pv,method='fdr')
    #RESULT$pcuse=PCUSE
    print('######################################')
    print('All Main Steps Finished !!!')
    print('######################################')
    return(RESULT)
    }



