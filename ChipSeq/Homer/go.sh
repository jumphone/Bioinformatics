#
#http://homer.ucsd.edu/homer/

perl configureHomer.pl -install homer
perl /Users/zha8dh/Desktop/CCHMC_Project/LIJUN/Motif_20190312/.//configureHomer.pl -install mm10

export PATH=$PATH:/Users/zha8dh/Desktop/CCHMC_Project/LIJUN/Motif_20190312/bin

findMotifsGenome.pl Lijun.txt.bed.bed mm10 OUTPUT -size 200

nohup findMotifsGenome.pl Lijun.txt.bed.POS.bed mm10 OUTPUT_POS -size 200 &

nohup findMotifsGenome.pl Lijun.txt.bed.NEG.bed mm10 OUTPUT_NEG -size 200 &

