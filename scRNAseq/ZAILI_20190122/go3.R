library(Seurat)
library(dplyr)

meta.data=read.table('cerebellum_cell_metadata.tsv',sep='\t',header=T)


TAG='E10A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E10A.data=as.matrix(exp.data[,used])

TAG='E10B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E10B.data=as.matrix(exp.data[,used])

TAG='E11A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E11A.data=as.matrix(exp.data[,used])

TAG='E11B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E11B.data=as.matrix(exp.data[,used])

TAG='E12A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E12A.data=as.matrix(exp.data[,used])

TAG='E12B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E12B.data=as.matrix(exp.data[,used])

TAG='E13A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E13A.data=as.matrix(exp.data[,used])

TAG='E13B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E13B.data=as.matrix(exp.data[,used])

TAG='E14A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E14A.data=as.matrix(exp.data[,used])

TAG='E14B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E14B.data=as.matrix(exp.data[,used])

TAG='E15A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E15A.data=as.matrix(exp.data[,used])

TAG='E15B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E15B.data=as.matrix(exp.data[,used])

TAG='E16A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E16A.data=as.matrix(exp.data[,used])

TAG='E16B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E16B.data=as.matrix(exp.data[,used])

TAG='E17A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E17A.data=as.matrix(exp.data[,used])

TAG='E17B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
E17B.data=as.matrix(exp.data[,used])

TAG='P00A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P00A.data=as.matrix(exp.data[,used])

TAG='P00B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P00B.data=as.matrix(exp.data[,used])

TAG='P04A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P04A.data=as.matrix(exp.data[,used])

TAG='P04B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P04B.data=as.matrix(exp.data[,used])

TAG='P07A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P07A.data=as.matrix(exp.data[,used])

TAG='P07B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P07B.data=as.matrix(exp.data[,used])

TAG='P10A'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P10A.data=as.matrix(exp.data[,used])

TAG='P10B'
exp.data <- Read10X(data.dir = paste0("./cerebellum_filtered_matrices/",TAG))
CN=colnames(exp.data)
NCN=paste0(TAG,'_',CN,'-1')
colnames(exp.data)=NCN
used=which(colnames(exp.data) %in% meta.data[,3])
P10B.data=as.matrix(exp.data[,used])



source('scRef.R')
E1=.simple_combine(E10A.data,E10B.data)
E2=.simple_combine(E11A.data,E11B.data)
E3=.simple_combine(E12A.data,E12B.data)
E4=.simple_combine(E13A.data,E13B.data)
E5=.simple_combine(E14A.data,E14B.data)
E6=.simple_combine(E15A.data,E15B.data)
E7=.simple_combine(E16A.data,E16B.data)
E8=.simple_combine(E17A.data,E17B.data)
E9=.simple_combine(P00A.data,P00B.data)
E10=.simple_combine(P04A.data,P04B.data)
E11=.simple_combine(P07A.data,P07B.data)
E12=.simple_combine(P10A.data,P10B.data)

rm(E10A.data)
rm(E10B.data)
rm(E11A.data)
rm(E11B.data)
rm(E12A.data)
rm(E12B.data)
rm(E13A.data)
rm(E13B.data)
rm(E14A.data)
rm(E14B.data)
rm(E15A.data)
rm(E15B.data)
rm(E16A.data)
rm(E16B.data)
rm(E17A.data)
rm(E17B.data)
rm(P00A.data)
rm(P00B.data)
rm(P04A.data)
rm(P04B.data)
rm(P07A.data)
rm(P07B.data)
rm(P10A.data)
rm(P10B.data)
rm(exp.data)

A1=.simple_combine(E1$combine,E2$combine)$combine
A2=.simple_combine(E3$combine,E4$combine)$combine
A3=.simple_combine(E5$combine,E6$combine)$combine
A4=.simple_combine(E7$combine,E8$combine)$combine
A5=.simple_combine(E9$combine,E10$combine)$combine
A6=.simple_combine(E11$combine,E12$combine)$combine


rm(E1)
rm(E2)
rm(E3)
rm(E4)
rm(E5)
rm(E6)
rm(E7)
rm(E8)
rm(E9)
rm(E10)
rm(E11)
rm(E12)


B1=.simple_combine(A1,A2)$combine
B2=.simple_combine(A3,A4)$combine
B3=.simple_combine(A5,A6)$combine


rm(A1)
rm(A2)
rm(A3)
rm(A4)
rm(A5)
rm(A6)


C1=.simple_combine(B1,B2)$combine
rm(B1)
rm(B2)

EXP=.simple_combine(C1,B3)$combine

rm(C1)
rm(B3)



CN=colnames(EXP)
MCN=meta.data[,3]
O=c()
i=1
for(one in MCN){
   O=c(O,which(CN ==one ))
   print(i);i=i+1
}

EXP=EXP[,O]

saveRDS(EXP,file='cerebullum_dev.RDS')










