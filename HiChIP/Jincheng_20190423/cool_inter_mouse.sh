
mkdir $1.tmp

cooler dump --join -r chr1 $1 > $1.tmp/chr1.txt
cooler dump --join -r chr2 $1 > $1.tmp/chr2.txt
cooler dump --join -r chr3 $1 > $1.tmp/chr3.txt
cooler dump --join -r chr4 $1 > $1.tmp/chr4.txt
cooler dump --join -r chr5 $1 > $1.tmp/chr5.txt
cooler dump --join -r chr6 $1 > $1.tmp/chr6.txt
cooler dump --join -r chr7 $1 > $1.tmp/chr7.txt
cooler dump --join -r chr8 $1 > $1.tmp/chr8.txt
cooler dump --join -r chr9 $1 > $1.tmp/chr9.txt
cooler dump --join -r chr10 $1 > $1.tmp/chr10.txt
cooler dump --join -r chr11 $1 > $1.tmp/chr11.txt
cooler dump --join -r chr12 $1 > $1.tmp/chr12.txt
cooler dump --join -r chr13 $1 > $1.tmp/chr13.txt
cooler dump --join -r chr14 $1 > $1.tmp/chr14.txt
cooler dump --join -r chr15 $1 > $1.tmp/chr15.txt
cooler dump --join -r chr16 $1 > $1.tmp/chr16.txt
cooler dump --join -r chr17 $1 > $1.tmp/chr17.txt
cooler dump --join -r chr18 $1 > $1.tmp/chr18.txt
cooler dump --join -r chr19 $1 > $1.tmp/chr19.txt
cooler dump --join -r chrX $1 > $1.tmp/chrX.txt
cooler dump --join -r chrY $1 > $1.tmp/chrY.txt
cooler dump --join -r chrM $1 > $1.tmp/chrM.txt


cat $1.tmp/chr* > $1.inter.txt
