library(MuSiC)
library('Biobase')
#GSE50244.bulk.eset = readRDS('GSE50244bulkeset.rds')


assayData=read.table('gene_rep_matix_anno_diff.txt.pure',sep='\t',row.names=1,header=T)
assayData=as.matrix(assayData)

sampleID=1:ncol(assayData)
SubjectName=colnames(assayData)
cellTypeID=colnames(assayData)
cellType=colnames(assayData)

pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
rownames(pheno.matrix)=colnames(assayData)
metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
Bulk.eset = ExpressionSet(assayData = data.matrix(assayData), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )


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
PP = music_prop(bulk.eset = Bulk.eset, sc.eset = SC.eset , clusters = 'cellType', samples = 'cellType',verbose = F)


