
mkdir $1.tmp

cooler dump --join -r chr1 $1 > $1.tmp/chr1.cool
cooler dump --join -r chr2 $1 > $1.tmp/chr2.cool
cooler dump --join -r chr3 $1 > $1.tmp/chr3.cool
cooler dump --join -r chr4 $1 > $1.tmp/chr4.cool
cooler dump --join -r chr5 $1 > $1.tmp/chr5.cool
cooler dump --join -r chr6 $1 > $1.tmp/chr6.cool
cooler dump --join -r chr7 $1 > $1.tmp/chr7.cool
cooler dump --join -r chr8 $1 > $1.tmp/chr8.cool
cooler dump --join -r chr9 $1 > $1.tmp/chr9.cool
cooler dump --join -r chr10 $1 > $1.tmp/chr10.cool
cooler dump --join -r chr11 $1 > $1.tmp/chr11.cool
cooler dump --join -r chr12 $1 > $1.tmp/chr12.cool
cooler dump --join -r chr13 $1 > $1.tmp/chr13.cool
cooler dump --join -r chr14 $1 > $1.tmp/chr14.cool
cooler dump --join -r chr15 $1 > $1.tmp/chr15.cool
cooler dump --join -r chr16 $1 > $1.tmp/chr16.cool
cooler dump --join -r chr17 $1 > $1.tmp/chr17.cool
cooler dump --join -r chr18 $1 > $1.tmp/chr18.cool
cooler dump --join -r chr19 $1 > $1.tmp/chr19.cool
cooler dump --join -r chrX $1 > $1.tmp/chrX.cool
cooler dump --join -r chrY $1 > $1.tmp/chrY.cool
cooler dump --join -r chrM $1 > $1.tmp/chrM.cool


cooler merge $1.inter.cool  $1.tmp/chr1.cool $1.tmp/chr2.cool $1.tmp/chr3.cool $1.tmp/chr4.cool $1.tmp/chr5.cool $1.tmp/chr6.cool  $1.tmp/chr7.cool $1.tmp/chr8.cool $1.tmp/chr9.cool $1.tmp/chr10.cool $1.tmp/chr11.cool $1.tmp/chr12.cool $1.tmp/chr13.cool $1.tmp/chr14.cool $1.tmp/chr15.cool $1.tmp/chr16.cool $1.tmp/chr17.cool $1.tmp/chr18.cool $1.tmp/chr19.cool $1.tmp/chrX.cool $1.tmp/chrY.cool $1.tmp/chrM.cool
