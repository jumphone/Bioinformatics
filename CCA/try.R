
source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/scPA.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
CT=readRDS('CT.RDS')
MS=readRDS('MS.RDS')
CTX=readRDS('CTX.RDS')
MSX=readRDS('MSX.RDS')

CTG=.getGroup(CTX,'CT',CNUM=100)
MSG=.getGroup(MSX,'MS',CNUM=100)

VP=.getValidpair(CT, CTG, MS, MSG, CPU=4, method='kendall', do.plot=FALSE, print_step=10)
VP=VP$vp


EXP=.simple_combine(CT,MS)$combine
GROUP=c(CTG,MSG)
CONDITION=c(rep('CT',ncol(CT)),rep('MS',ncol(MS)))


library(Seurat)
pbmc=CreateSeuratObject(raw.data = EXP, min.cells = 0, min.genes = 0, project = "ALL")
pbmc@meta.data$group=GROUP
pbmc@meta.data$condition=CONDITION

MAP=rep('NA',length(GROUP))
MAP[which(GROUP %in% VP[,1])]='CT'
MAP[which(GROUP %in% VP[,2])]='MS'
pbmc@meta.data$map=MAP

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"))


PCNUM=50
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)

saveRDS(pbmc, file='try_pbmc.RDS')


DR=pbmc@dr$pca@cell.embeddings
saveRDS(DR, file='try_DR.RDS')


B1index=which(CONDITION=='CT')
B2index=which(CONDITION=='MS')


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
            if(x1 <=min_lst1){x1=min_lst1}
            if(x1 >=max_lst1){x1=max_lst1}
                
            x1=x1
            dlst1=c()
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst1[i]
                this_mean=mean_lst1[i]
                this_d=dnorm(x1,sd=this_sd,mean=this_mean)
                dlst1=c(dlst1,this_d)
                i=i+1} 
                out=sum(dlst1/sum(dlst1)*mean_com)
            return(out)}
      
        .x2_to_com=function(x2){
            if(x2 <=min_lst2){x2=min_lst2}
            if(x2 >=max_lst2){x2=max_lst2}
                
            x2=x2
            dlst2=c()
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst2[i]
                this_mean=mean_lst2[i]
                this_d=dnorm(x2,sd=this_sd,mean=this_mean)
                dlst2=c(dlst2,this_d)
                i=i+1} 
                out=sum(dlst2/sum(dlst2)*mean_com)
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



#approxfun(density(THIS_PC[p1]




DR=pbmc@dr$pca@cell.embeddings
B1index=which(CONDITION=='CT')
B2index=which(CONDITION=='MS')

OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)








DR=pbmc@dr$pca@cell.embeddings
B1index=which(CONDITION=='CT')
B2index=which(CONDITION=='MS')

OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)


pbmc@dr$oldpca=pbmc@dr$pca
pbmc@dr$pca@cell.embeddings=OUT$adr


PCUSE=which(p.adjust(OUT$pv,method='fdr')<0.05 & OUT$cor>0.8)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE)

