KO=P5-50K-jincheng.rmdup.bam
WT=PDGFRa-GFP-50K-E14.5-jincheng.rmdup.bam

cat *_peaks.* > ALL_PEAK.bed
bedtools sort -i ALL_PEAK.bed > ALL_PEAK.sorted.bed
bedtools merge -i ALL_PEAK.sorted.bed > ALL_PEAK.sorted.merged.bed
samtools stats $KO > $KO\.stats
samtools stats $WT > $WT\.stats 
python addname.py
bedtools coverage -abam $KO -b ALL_PEAK.sorted.merged.named.bed > $KO\.ALLPEAK 
bedtools coverage -abam $WT -b ALL_PEAK.sorted.merged.named.bed > $WT\.ALLPEAK 
python combine_peak_RPKM.py
python rm20.py

