setwd('F:/Zhenglab/Combine')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')


REF=.readTable(PATH='GSE92332_atlas_UMIcounts.txt',SEP='\t')


getbatch <- function(x){
    y=unlist(strsplit(x, "_"))
    y=y[length(y)]
    return(y)
}
CN=colnames(REF)
BATCH=apply(matrix(CN,ncol=1),1,getbatch)
table(BATCH)
LABEL=BATCH

library(dplyr)
library(Seurat)

#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

CDC42HET <- Read10X(data.dir = "./CDC42_HET")
CDC42KO<- Read10X(data.dir = "./Small_Intestine_KO")
AGE <- Read10X(data.dir = "./age")
YOUNG <- Read10X(data.dir = "./young")

BATCH=c(rep('NATURE',ncol(REF)),
        rep('CDC42HET',ncol(CDC42HET)),
        rep('CDC42KO',ncol(CDC42KO)),   
        rep('AGE',ncol(AGE)),
       rep('YOUNG',ncol(YOUNG))
       )
####################

D1=.simple_combine(REF, CDC42HET)$combine
rm(REF)
rm(CDC42HET)
gc()


D2=.simple_combine(CDC42KO, AGE)$combine
rm(CDC42KO)
rm(AGE)
gc()


D3=.simple_combine(D1, D2)$combine
rm(D1)
rm(D2)
gc()

DATA=.simple_combine(D3, YOUNG)$combine
rm(D3)
rm(YOUNG)
gc()
saveRDS(DATA,'DATA.RDS')
saveRDS(BATCH,'BATCH.RDS')


TMP=rep('NA',ncol(DATA))
TMP[which(BATCH=='NATURE')]=LABEL
LABEL=TMP
saveRDS(BATCH,'LABEL.RDS')


###################################################################

getNon0=function(x){
   return(length(which(x>0)))
   }

NON0=apply(DATA,2,getNon0)

NON0.G=apply(DATA,1,getNon0)







