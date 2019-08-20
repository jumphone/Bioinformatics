

setwd('/users/zha8dh/tianlab/HFZ')
source('./BEER.R')
DATA=readRDS('./DATA.RDS')
BATCH=readRDS('./BATCH.RDS')

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   
saveRDS(mybeer,file='mybeer.RDS')


##################################



setwd('/users/zha8dh/tianlab/HFZ')
source('./BEER.R')
.set_python('')
mybeer=readRDS('./mybeer.RDS')





