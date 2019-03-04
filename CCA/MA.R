
#library(Seurat)
#source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


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


.getValidpair <- function(DATA1, GROUP1, DATA2, GROUP2, CUTOFF=0, do.plot=TRUE){
    source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('Start')
    print('Step1.Generate Reference...')
    REF1=.generate_ref(DATA1, cbind(GROUP1, GROUP1), min_cell=1) 
    REF2=.generate_ref(DATA2, cbind(GROUP2, GROUP2), min_cell=1) 
    print('Step2.Calculate Correlation Coefficient...')
    out = .get_cor( REF1, REF2, method='kendall',CPU=4, print_step=10)
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
    if(do.plot==TRUE){plot(C)}
    VP=VP[which(C>=CUTOFF),]  
    print('Finished!!!')
    return(VP)
    }










