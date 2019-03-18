load('GliomaSingleCell.RData')


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

mybeer=MBEER(DATA, BATCH, MAXBATCH="", CNUM=10, PCNUM=20, CPU=2, SEED=1 )

