
library(Seurat)
source('scRef.R')

tag_data=read.table('GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt',sep='\t',row.names=1,header=T)


pbmc.data=read.table('GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt',sep='\t',row.names=1,header=T)
#saveRDS(pbmc.data,'ALL.RDS')

MS=pbmc.data[,which(tag_data[,3]=='MS')]
CT=pbmc.data[,which(tag_data[,3]=='Ctrl')]

saveRDS(MS,'MS.RDS')
saveRDS(CT,'CT.RDS')

CPU=3
PCNUM=50
PCUSE=1:PCNUM

##########
MSD = CreateSeuratObject(raw.data = MS, min.cells = 0, min.genes = 0, project = "MS") 
MSD <- NormalizeData(object = MSD, normalization.method = "LogNormalize", scale.factor = 10000)
MSD <- ScaleData(object = MSD, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
MSD <- RunPCA(object = MSD, pcs.compute=PCNUM, pc.genes = rownames(MSD@data), do.print = FALSE)
MSD <- RunTSNE(object = MSD, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
MSX=MSD@dr$tsne@cell.embeddings
saveRDS(MSX,file='MSX.RDS')
##########
CTD = CreateSeuratObject(raw.data = CT, min.cells = 0, min.genes = 0, project = "CT") 
CTD <- NormalizeData(object = CTD, normalization.method = "LogNormalize", scale.factor = 10000)
CTD <- ScaleData(object = CTD, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
CTD <- RunPCA(object = CTD, pcs.compute=PCNUM, pc.genes = rownames(CTD@data), do.print = FALSE)
CTD <- RunTSNE(object = CTD, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
CTX=CTD@dr$tsne@cell.embeddings
saveRDS(CTX,file='CTX.RDS')


##########
##########
##########
##########
##########
##########
library(Seurat)
source('scRef.R')

MS=readRDS("MS.RDS")
MSX=readRDS("MSX.RDS")
CT=readRDS("CT.RDS")
CTX=readRDS("CTX.RDS")




.getGroup=function(X,GNUM){
    GNUM=GNUM
    DR=X
    RANK=rank(DR,ties.method='random')
    CUTOFF=round(max(RANK)/GNUM)
    GROUP=rep('NA',length(RANK))
    i=1
    j=1
    while(i<=length(MSR)){
        GROUP[which(RANK==i)]=paste0('G_',as.character(j))
        if(i%%CUTOFF==1){j=j+1}
        i=i+1}
    return(GROUP)
}



#####
GNUM=100
MSG=.getGroup(MSX,GNUM)
CTG=.getGroup(CTX,GNUM)
########
MSR=.generate_ref(MS, cbind(MSG,MSG), min_cell=1) 
CTR=.generate_ref(CT, cbind(CTG,CTG), min_cell=1) 

saveRDS(MSR,file='MSR.RDS')
saveRDS(CTR,file='CTR.RDS')
########

out = .get_cor( MSR, CTR, method='kendall',CPU=4, print_step=10)











