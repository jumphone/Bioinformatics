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
