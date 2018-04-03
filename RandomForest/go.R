#install.packages("randomForest")

library(randomForest)

rawdata=read.table('OMIMID_Name.txt.stat.name.sig.t',header=T,row.names=1)
trawdata=t(rawdata)

TAG=read.table('OMIMID_Name.txt.stat.name.sig.t.tag',header=FALSE)[,1]
TARGET=as.character(TAG) #as.factor(TAG)


#DATA=data.frame(trawdata)
#DATA$TARGET=TARGET

set.seed(12345)
index <- sample(2,nrow(DATA),replace = TRUE,prob=c(0.7,0.3))
traindata <- trawdata[index==1,]
testdata <- trawdata[index==2,]

traintarget <- TARGET[index==1]
traintarget=as.factor(traintarget)

testtarget <- TARGET[index==2]

set.seed(12345)
rf_ntree <- randomForest(traintarget~. ,data=traindata,ntree=500,proximity=TRUE)

pdf('tmp.pdf')
plot(rf_ntree)
dev.off()

rf_ntree_pred <- predict(rf_ntree, newdata=testdata)
table(rf_ntree_pred, testtarget)


