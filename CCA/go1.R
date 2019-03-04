source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/CCA/scPA.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

CT=readRDS('CT.RDS')
MS=readRDS('MS.RDS')

CTX=.data2one(CT)
MSX=.data2one(MS)


saveRDS(CT,file='CT.RDS')
saveRDS(MS,file='MS.RDS')
saveRDS(CTX,file='CTX.RDS')
saveRDS(MSX,file='MSX.RDS')


CTG=.getGroup(CTX,'CT',CNUM=100)
MSG=.getGroup(MSX,'MS',CNUM=100)

VP=.getValidpair(CT, CTG, MS, MSG, CPU=4, method='kendall', do.plot=FALSE, print_step=10)
  

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
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"))

PCNUM=100
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)


DR=pbmc@dr$pca@cell.embeddings
B1index=which(CONDITION=='CT')
B2index=which(CONDITION=='MS')



.dr2adr <- function(DR, B1index, B2index, GROUP, VP){
    library(dtw)
    library(MALDIquant)
    library(pcaPP)
    OUT=list()
    OUT$adr=DR
    VALID_PAIR=VP
    ALL_COR=c()   
    ALL_PW=c() 
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
        comlst1o=order(comlst[,1])
        comlst2o=order(comlst[,2])
        
        lst1lst1=compc[comlst1o][match.closest(DR[index1,THIS_DR], comlst[,1][comlst1o])]
        lst2lst2=compc[comlst2o][match.closest(DR[index2,THIS_DR], comlst[,2][comlst2o])]
        
        OUT$adr[index1,THIS_DR]=lst1lst1
        OUT$adr[index2,THIS_DR]=lst2lst2
        

        this_pw=var(DR[c(vindex1,vindex2),THIS_DR])#/sd(DR[,THIS_DR])
        ALL_PW=c(ALL_PW, this_pw)
      
        #this_cor=cor.fk(OUT$adr[,THIS_DR],DR[,THIS_DR])
        #ALL_COR=c(ALL_COR, this_cor)
        
        #lst1tocom = function(x){quantile(compc,ecdf(comlst[,1])(x))}   
        #lst2tocom = function(x){quantile(compc,ecdf(comlst[,2])(x))}
         
        #OUT$adr[index1,THIS_DR]=lst1tocom(DR[index1,THIS_DR])
        #OUT$adr[index2,THIS_DR]=lst2tocom(DR[index2,THIS_DR])
         
        print(THIS_DR)
        THIS_DR=THIS_DR+1}
    
    #OUT$cor=ALL_COR
    OUT$var=ALL_PW
    print('Finished!!!')
    return(OUT)
    }



OUT=.dr2adr(DR, B1index, B2index, GROUP, VP)
plot(OUT$var,type='l')




pbmc@dr$alnpca=pbmc@dr$pca
pbmc@dr$alnpca@key='APC'
pbmc@dr$alnpca@cell.embeddings=OUT$adr


PCUSE=which(OUT$var>mean(OUT$var))
pbmc <- RunTSNE(object = pbmc, reduction.use='alnpca',dims.use = PCUSE, do.fast = TRUE)

DimPlot(object =pbmc, reduction.use = "tsne", group.by = "map",  pt.size = 0.5, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "tsne", group.by = "condition",  pt.size = 0.5, do.return = TRUE)
DimPlot(object =pbmc, reduction.use = "tsne",  pt.size = 0.5, do.return = TRUE)






