setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')


REF=.readTable(PATH='GSE92332_AtlasFullLength_TPM.txt',SEP='\t')


getbatch <- function(x){
    y=unlist(strsplit(x, "_"))
    y=y[length(y)]
    return(y)
}
CN=colnames(REF)
BATCH=apply(matrix(CN,ncol=1),1,getbatch)
