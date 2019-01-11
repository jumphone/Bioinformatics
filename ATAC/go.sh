KO=
WT=

/usr/bin/macs2 callpeak -t  $WT $KO -n ./POOL -B -f BAM --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01 -g mm

cat *_peaks.narrowPeak > ALL_PEAK.bed
bedtools sort -i ALL_PEAK.bed > ALL_PEAK.sorted.bed
bedtools merge -i ALL_PEAK.sorted.bed > ALL_PEAK.sorted.merged.bed
samtools stats $KO > $KO\.stats
samtools stats $WT > $WT\.stats
python addname.py
bedtools coverage -abam $KO -b ALL_PEAK.sorted.merged.named.bed > $KO\.ALLPEAK
bedtools coverage -abam $WT -b ALL_PEAK.sorted.merged.named.bed > $WT\.ALLPEAK
python combine_peak_RPKM.py $KO\.ALLPEAK $WT\.ALLPEAK
python rm20.py
python getgene.py
