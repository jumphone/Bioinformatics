#install.packages("randomForest")

library(randomForest)

rawdata=read.table('OMIMID_Name.txt.stat.name.sig.t',header=T,row.names=1)
trawdata=t(rawdata)

TAG=read.table('OMIMID_Name.txt.stat.name.sig.t.tag',header=F)[,1]
TARGET_STR=TAG
TARGET_NUM=as.numeric(as.factor(TAG))

DATA=cbind(trawdata,TARGET_NUM)

set.seed(12345)
index <- sample(2,nrow(DATA),replace = TRUE,prob=c(0.7,0.3))
traindata <- DATA[index==1,]
testdata <- DATA[index==2,]

set.seed(12345)
rf_ntree <- randomForest(TARGET_NUM~.,data=traindata,ntree=100,proximity=TRUE)

pdf('tmp.pdf')
plot(rf_ntree)
dev.off()

rf_ntree_pred <- predict(rf_ntree, newdata=testdata)
table(rf_ntree_pred, testdata$TARGET_NUM)


