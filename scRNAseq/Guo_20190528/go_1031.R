setwd('F:/GuoFuKun/NEW1031/GSE114727_RAW/')



###################################
PATH="GSM3148640_BC08_TUMOR3_counts.csv.gz"
TAG='TUMOR3.BC08'
################################
this_data=read.csv(gzfile(PATH), sep=',',header=T)
this_data=this_data[,c(2:ncol(this_data))]
this_data=t(this_data)
this_data[which(is.na(this_data))]=0

colnames(this_data)=paste0(TAG,'_',c(1:ncol(this_data)))
this_data[1:3,1:3]
dim(this_data)
###############################
D1=this_data
#######################



