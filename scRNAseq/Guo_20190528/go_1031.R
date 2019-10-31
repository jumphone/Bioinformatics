setwd('F:/GuoFuKun/NEW1031/GSE114727_RAW/')


.readData=function(PATH,TAG){
    PATH=PATH
    TAG=TAG
    this_data=read.csv(gzfile(PATH), sep=',',header=T)
    this_data=this_data[,c(2:ncol(this_data))]
    this_data=t(this_data)
    this_data[which(is.na(this_data))]=0

    colnames(this_data)=paste0(TAG,'_',c(1:ncol(this_data)))
    print(this_data[1:3,1:3])
    print(dim(this_data))
    return(this_data)

    }


###################################
PATH="GSM3148640_BC08_TUMOR3_counts.csv.gz"
TAG='TUMOR3.BC08'
################################
D1=.readData(PATH,TAG)
#######################


###################################
PATH="GSM3148639_BC08_TUMOR2_counts.csv.gz"
TAG='TUMOR2.BC08'
################################
this_data=read.csv(gzfile(PATH), sep=',',header=T)
this_data=this_data[,c(2:ncol(this_data))]
this_data=t(this_data)
this_data[which(is.na(this_data))]=0

colnames(this_data)=paste0(TAG,'_',c(1:ncol(this_data)))
this_data[1:3,1:3]
dim(this_data)
###############################
D2=this_data
#######################


###################################
PATH="GSM3148638_BC08_TUMOR1_counts.csv.gz"
TAG='TUMOR2.BC08'
################################
this_data=read.csv(gzfile(PATH), sep=',',header=T)
this_data=this_data[,c(2:ncol(this_data))]
this_data=t(this_data)
this_data[which(is.na(this_data))]=0

colnames(this_data)=paste0(TAG,'_',c(1:ncol(this_data)))
this_data[1:3,1:3]
dim(this_data)
###############################
D3=this_data
#######################













