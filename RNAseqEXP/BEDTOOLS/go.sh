
GENE=./HumanHg19GeneSort.bed
BAM=XXX.bam
sam2bed="python2.7 ./sam2bed.py"
cov2rpkm="python2.7 ./cov2rpkm.py"

samtools view $BAM > $BAM\.sam
$sam2bed $BAM\.sam $BAM\.sam.bed
bedtools coverage -a $GENE -b $BAM\.sam.bed -wa -wb > $BAM\.gene.cov
$cov2rpkm $BAM\.gene.cov $BAM\.gene.rpkm.txt
