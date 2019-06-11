library(MuSiC)
library('Biobase')
#GSE50244.bulk.eset = readRDS('GSE50244bulkeset.rds')


assayData=read.table('gene_rep_matix_anno_diff.txt.pure',sep='\t',row.names=1,header=T)
assayData=as.matrix(assayData)

#sampleID=1:ncol(assayData)
#SubjectName=colnames(assayData)
#cellTypeID=colnames(assayData)
#cellType=colnames(assayData)

#sampleType=colnames(assayData)
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



REF=.generate_ref(COM,cbind(CELLTYPE,CELLTYPE),M='MEAN')
VAR=apply(REF,1,var)
REF=REF[which(VAR> quantile(VAR,0.9)),]

write.table(REF,file='cibersort_sig.txt',sep='\t',row.names=T,col.names=T,quote=F)





sampleID=1:ncol(COM)
SubjectName=colnames(COM)
cellTypeID=CELLTYPE
cellType=CELLTYPE
pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
rownames(pheno.matrix)=colnames(COM)
metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
SC.eset = ExpressionSet(assayData = data.matrix(COM), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )
#SC.eset=ExpressionSet(COM)

#sampleType=c(rep('E14',ncol(D1)),rep('P0',ncol(D2)))
sampleType=colnames(COM)#rep('ALL',ncol(COM))

library(xbioc)
library(pvar)
PP_ALL = music_prop(bulk.eset = Bulk.eset, sc.eset = SC.eset ,markers = NULL, clusters = cellType, samples = sampleType, verbose = F)



saveRDS(PP_ALL,file='PP_ALL.RDS')






OUT=PP_ALL$Est.prop.weighted
S=apply(OUT,2,sum)
OUT=OUT[,which(S>0)]




library('gplots')

pdf('HEAT.pdf',width=7,height=5)
heatmap.2(OUT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F, trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(15,5))
heatmap.2(OUT,scale=c("column"),dendrogram='none',Colv=F,Rowv=F, trace='none',col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(15,5))
dev.off()


write.table(OUT, file='MUSIC_OUT.txt',sep='\t',quote=F,row.names=T,col.names=T)


















































