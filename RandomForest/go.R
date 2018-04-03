#install.packages("randomForest")

library(randomForest)

rawdata=read.table('OMIMID_Name.txt.stat.name.sig.t',header=T,row.names=1)
trawdata=t(rawdata)

TAG=read.table('OMIMID_Name.txt.stat.name.sig.t.tag',header=T)[,1]

set.seed(12345)
index <- sample(2,nrow(trawdata),replace = TRUE,prob=c(0.7,0.3))
traindata <- trawdata[index==1,]
testdata <- trawdata[index==2,]

TARGET=rep(0,length(traindata[,1]))
TARGET[which(TAG=='OMIM220210')]=1

set.seed(12345)
rf_ntree <- randomForest(TARGET~.,data=traindata,ntree=300)

pdf('tmp.pdf')
plot(rf_ntree)
dev.off()
