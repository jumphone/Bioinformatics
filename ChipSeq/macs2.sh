bam=$1
species=$2 #mm,hs
macs2 callpeak -t $bam -f BAM -g $species -n $bam\.macs2_narrow -B -q 0.01
macs2 callpeak -t $bam  --broad -g $species -n $bam\.macs2_broad --broad-cutoff 0.1
