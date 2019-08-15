
#setwd('F:/Zhenglab/Combine')

setwd('/Volumes/Feng/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#REF=.readTable(PATH='GSE92332_atlas_UMIcounts.txt',SEP='\t')
REF=readRDS('REF.RDS')
getbatch <- function(x){
    y=unlist(strsplit(x, "_"))
    y=y[length(y)]
    return(y)
}

#saveRDS(REF,file='REF.RDS')


CN=colnames(REF)
BATCH=apply(matrix(CN,ncol=1),1,getbatch)
table(BATCH)
LABEL=BATCH




############################

CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO<- Read10X(data.dir = "./Small_Intestine_KO")
AGE <- Read10X(data.dir = "./age")
YOUNG <- Read10X(data.dir = "./young")

BATCH1=c(rep('NATURE',ncol(REF)),
        rep('CDC42HET',ncol(CDC42HET)),
        rep('CDC42KO',ncol(CDC42KO))
       )

BATCH2=c(rep('NATURE',ncol(REF)),
        rep('AGE',ncol(AGE)),
       rep('YOUNG',ncol(YOUNG))
)
####################

D1=.simple_combine(REF, CDC42HET)$combine
rm(CDC42HET)
gc()


DATA1=.simple_combine(D1, CDC42KO)$combine
rm(CDC42KO)
gc()

D2=.simple_combine(REF, AGE)$combine
rm(AGE)
gc()


DATA2=.simple_combine(D2, YOUNG)$combine
rm(YOUNG)
gc()

####################
saveRDS(DATA1,'DATA1.RDS') # Zheng Zhang
saveRDS(BATCH1,'BATCH1.RDS')

saveRDS(DATA1,'DATA2.RDS') # Age Young
saveRDS(BATCH1,'BATCH2.RDS')
####################















