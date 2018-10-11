F1=$1
F2=$2
O=$3

GENE=./Mus_musculus.GRCm38.87.chr.gtf.combined.pc.bed

mkdir $O

cat $F1 $F2 > $O\/combine.bed
bedtools sort -i $O\/combine.bed > $O\/combine.sorted.bed
bedtools merge -i $O\/combine.sorted.bed > $O\/combine.sorted.merged.bed
python gettag.py $O\/combine.sorted.merged.bed $O\/combine.sorted.merged.taged.bed
bedtools coverage -a $O\/combine.sorted.merged.taged.bed -b $F1 > $O\/combine.sorted.merged.taged.bed.F1.cov
bedtools coverage -a $O\/combine.sorted.merged.taged.bed -b $F2 > $O\/combine.sorted.merged.taged.bed.F2.cov
python combine.py $O\/combine.sorted.merged.taged.bed $O\/combine.sorted.merged.taged.bed.F1.cov $O\/combine.sorted.merged.taged.bed.F2.cov $O\/compare.txt
bedtools intersect -a $O\/compare.txt -b $GENE -wa -wb > $O\/compare.gene.txt
Rscript drawVENN.R $O 
