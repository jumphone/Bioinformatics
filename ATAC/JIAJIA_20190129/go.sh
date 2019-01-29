WT=./mOPC_ctrl_RL1.rmdup.bam
KO=./mOPC_Brg1cKO_RL2.rmdup.bam


macs2 callpeak -t  $WT $KO -n ./POOL -B -f BAM --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01 -g mm
cat *_peaks.narrowPeak > ALL_PEAK.bed
bedtools sort -i ALL_PEAK.bed > ALL_PEAK.sorted.bed
bedtools merge -i ALL_PEAK.sorted.bed > ALL_PEAK.sorted.merged.bed
python addname.py
bedtools coverage -abam $KO -b ALL_PEAK.sorted.merged.named.bed > $KO\.ALLPEAK 
bedtools coverage -abam $WT -b ALL_PEAK.sorted.merged.named.bed > $WT\.ALLPEAK 

python rmNA.py $KO\.ALLPEAK
python rmNA.py $WT\.ALLPEAK

bedtools intersect -wa -wb -b $KO\.ALLPEAK.rmNa -a Mus_musculus.GRCm38.87.chr.gtf.combined.pc.bed.sort  > $KO\.ALLPEAK.gene
bedtools intersect -wa -wb -b $WT\.ALLPEAK.rmNa -a Mus_musculus.GRCm38.87.chr.gtf.combined.pc.bed.sort  > $WT\.ALLPEAK.gene
python getRPKM.py $KO\.ALLPEAK.gene $KO\.ALLPEAK.gene.rpkm
python getRPKM.py $WT\.ALLPEAK.gene $WT\.ALLPEAK.gene.rpkm

Rscript combineR.R $WT\.ALLPEAK.gene.rpkm $KO\.ALLPEAK.gene.rpkm
