
GENE=/users/zha8dh/tianlab/GBMLGG/workspace/HumanHg19GeneSort.bed
BAM=/users/zha8dh/tianlab/GBMLGG/data/228acc60-d652-4f33-86ce-4851cc271edf/UNCID_1415217.106c0a27-64c4-4496-83f3-1c05df492fd2.sorted_genome_alignments.bam
sam2bed="python2.7 /users/zha8dh/tianlab/GBMLGG/workspace/sam2bed.py"
cov2rpkm="python2.7 /users/zha8dh/tianlab/GBMLGG/workspace/cov2rpkm.py"

samtools view $BAM > $BAM\.sam
$sam2bed $BAM\.sam $BAM\.sam.bed
bedtools coverage -a $GENE -b $BAM\.sam.bed -wa -wb > $BAM\.gene.cov
$cov2rpkm $BAM\.gene.cov $BAM\.gene.rpkm.txt
