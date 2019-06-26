cat *peaks* > ALL.peak.bed
bedtools sort -i ALL.peak.bed > ALL.peak.sort.bed
bedtools merge -i ALL.peak.sort.bed >  ALL.peak.sort.merge.bed

BAM=Lijun_ATAC_1.rmdup.bam
bedtools bamtobed -i $BAM > $BAM\.bed
bedtools coverage -a ALL.peak.sort.merge.bed -b $BAM\.bed  > $BAM\.bed.cov
python rpkm.py $BAM\.bed.cov


BAM=Lijun_ATAC_2.rmdup.bam
bedtools bamtobed -i $BAM > $BAM\.bed
bedtools coverage -a ALL.peak.sort.merge.bed -b $BAM\.bed  > $BAM\.bed.cov
python rpkm.py $BAM\.bed.cov

BAM=Lijun_ATAC_3.rmdup.bam
bedtools bamtobed -i $BAM > $BAM\.bed
bedtools coverage -a ALL.peak.sort.merge.bed -b $BAM\.bed  > $BAM\.bed.cov
python rpkm.py $BAM\.bed.cov

BAM=Lijun_ATAC_4.rmdup.bam
bedtools bamtobed -i $BAM > $BAM\.bed
bedtools coverage -a ALL.peak.sort.merge.bed -b $BAM\.bed  > $BAM\.bed.cov
python rpkm.py $BAM\.bed.cov


bedtools intersect -a ALL.peak.sort.merge.bed -b rn5_refseq.txt.bed -wa -wb > INTER.txt
bedtools intersect -a ALL.peak.sort.merge.bed -b rn5_refseq.txt.bed.5kb.bed -wa -wb > INTER5kb.txt


python combine.py
python combine5kb.py