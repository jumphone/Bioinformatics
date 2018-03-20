cat *_peaks.bed > ALL.peak
python /home/zhangfeng/tools/bed_sort.py ALL.peak
bedtools merge -i ALL.peak.sort > ALL.peak.sort.merged

WP=../mapping_mm10/

WT=$WP\LK-wt-myc.rmdup.bam #HSC-wt-myc-2.rmdup.bam #HSC-WT-K27Ac.rmdup.bam
grep 'reads mapped' $WT\.stats
KO=$WP\LK-KO-myc.rmdup.bam #HSC-KO-Myc-2.rmdup.bam #HSC-KO-K27Ac.rmdup.bam
grep 'reads mapped' $KO\.stats
KI=$WP\LK-KI-myc.rmdup.bam #HSC-KI-Myc-2.rmdup.bam #HSC-KI-K27Ac.rmdup.bam
grep 'reads mapped' $KI\.stats

nohup bedtools coverage -abam $WT -b ALL.peak.sort.merged > ALL.peak.sort.merged.WT.cov &
nohup bedtools coverage -abam $KO -b ALL.peak.sort.merged > ALL.peak.sort.merged.KO.cov &
nohup bedtools coverage -abam $KI -b ALL.peak.sort.merged > ALL.peak.sort.merged.KI.cov &
