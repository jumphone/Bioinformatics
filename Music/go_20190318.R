library(MuSiC)
library('Biobase')
#GSE50244.bulk.eset = readRDS('GSE50244bulkeset.rds')


assayData=read.table('gene_rep_matix_anno_diff.txt.pure',sep='\t',row.names=1,header=T)
assayData=as.matrix(assayData)

#sampleID=1:ncol(assayData)
#SubjectName=colnames(assayData)
#cellTypeID=colnames(assayData)
#cellType=colnames(assayData)

sampleType=colnames(assayData)
#sampleID=colnames(assayData)
Bulk.eset=ExpressionSet(assayData)

#pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
#rownames(pheno.matrix)=colnames(assayData)
#metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
#Bulk.eset = ExpressionSet(assayData = data.matrix(assayData), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )







#####################

D1=read.table('GSE123335_E14_combined_matrix.txt',sep='\t',row.names=1,header=T)
D2=read.table('GSE123335_P0_combined_matrix.txt',sep='\t',row.names=1,header=T)
saveRDS(D1,file='D1.RDS')
saveRDS(D2,file='D2.RDS')

T1=read.table('E14_combined_matrix_ClusterAnnotations.txt',sep='\t',row.names=1,header=T)
T1[,1]=as.character(T1[,1])
AT1=c()
i=1
while(i<=ncol(D1)){
this_t='NA'
if(colnames(D1)[i] %in% rownames(T1)){
    this_t=T1[which(rownames(T1)==colnames(D1)[i]),1]    
    }
AT1=c(AT1,this_t)
i=i+1
}


T2=read.table('P0_combined_matrix_ClusterAnnotations.txt',sep='\t',row.names=1,header=T)
T2[,1]=as.character(T2[,1])
AT2=c()
i=1
while(i<=ncol(D2)){
this_t='NA'
if(colnames(D2)[i] %in% rownames(T2)){
    this_t=T2[which(rownames(T2)==colnames(D2)[i]),1]    
    }
AT2=c(AT2,this_t)
i=i+1
}



source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
COM=.simple_combine(D1,D2)$combine
COM=as.matrix(COM)
CELLTYPE=c(AT1,AT2)


#sampleID=1:ncol(COM)
#SubjectName=colnames(COM)
#cellTypeID=CELLTYPE
cellType=CELLTYPE
#pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
#rownames(pheno.matrix)=colnames(COM)
#metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
#SC.eset = ExpressionSet(assayData = data.matrix(COM), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )
SC.eset=ExpressionSet(COM)

sampleType=c(rep('E14',ncol(D1)),rep('P0',ncol(D2)))


library(xbioc)
library(pvar)
PP_ALL = music_prop(bulk.eset = Bulk.eset, sc.eset = SC.eset ,markers = NULL, clusters = cellType, samples = sampleType, verbose = F)







































#####################
library(data.table)
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

d1=read.table('Brain_ref_mouse.txt',sep='\t',row.names=1,header=T)
d2=read.table('FetalBrain_ref_mouse.txt',sep='\t',row.names=1,header=T)
D=.simple_combine(d1,d2)$combine
exp_ref_mat=D
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
NewRef=as.matrix(NewRef)
oldheader=colnames(NewRef)
newheader=c(paste0(colnames(NewRef),'1'),paste0(colnames(NewRef),'2'))
NewRef=cbind(NewRef,jitter(NewRef))
colnames(NewRef)=newheader


sampleID=1:ncol(NewRef)
SubjectName=colnames(NewRef)
cellTypeID=colnames(NewRef)
cellType=c(oldheader,oldheader)
pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
rownames(pheno.matrix)=colnames(NewRef)
metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
SC.eset = ExpressionSet(assayData = data.matrix(NewRef), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )

library(xbioc)
library(pvar)
PP_ALL = music_prop(bulk.eset = Bulk.eset, sc.eset = SC.eset , clusters = cellType, samples = sampleID, verbose = F)





































#####################
library(data.table)
NewRef=read.table('Zeisel_exp_sc_mat.txt',sep='\t',row.names=1,header=T)
NewRef=as.matrix(NewRef)
TAG=read.table('Zeisel_TAG.txt',sep='\t',row.names=1,header=T)
TAG=as.character(TAG[,1])
sampleID=1:ncol(NewRef)
SubjectName=colnames(NewRef)
cellTypeID=TAG
cellType=TAG
pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
rownames(pheno.matrix)=colnames(NewRef)
metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
SC.eset = ExpressionSet(assayData = data.matrix(NewRef), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )

library(xbioc)
library(pvar)
PP = music_prop(bulk.eset = Bulk.eset, sc.eset = SC.eset , clusters = cellType, samples = sampleID, verbose = F)





