D1=readRDS(file='D1.RDS')
D2=readRDS(file='D2.RDS')

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
