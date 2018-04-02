#install.packages("randomForest")

library(randomForest)
rawdata=read.table('OMIMID_Name.txt.stat.name.sig.t',header=T,row.names=1)
trawdata=t(rawdata)
index <- sample(2,nrow(iris),replace = TRUE,prob=c(0.7,0.3))




