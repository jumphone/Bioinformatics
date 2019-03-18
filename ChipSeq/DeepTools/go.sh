#BAM (sorted) to BigWig (BW)

INPUT=CB_P30_CB_H3k27ac_In3.rmdup.bam
samtools index $INPUT
bamCoverage -b $INPUT -o $INPUT\.bw

INPUT=CB_P30_CB_Olig2_In2.rmdup.bam
samtools index $INPUT
bamCoverage -b $INPUT -o $INPUT\.bw

INPUT=MB_P15_MB_H3K27ac_In7.rmdup.bam
samtools index $INPUT
bamCoverage -b $INPUT -o $INPUT\.bw

INPUT=MB_Ptch_Olig2_2mm10.bam
samtools index $INPUT
bamCoverage -b $INPUT -o $INPUT\.bw



CB_H3K27AC=CB_P30_CB_H3k27ac_In3.rmdup.bam.bw
CB_OLIG2=CB_P30_CB_Olig2_In2.rmdup.bam.bw
MB_H3K27AC=MB_P15_MB_H3K27ac_In7.rmdup.bam.bw
MB_OLIG2=MB_Ptch_Olig2_2mm10.bam.bw

CB_H3K27AC_summit=./H3k27ac_CB/H3k27ac_CB_summits.bed
CB_OLIG2_summit=./olig2_CB/Olig2_CB_summits.bed
MB_H3K27AC_summit=./H3k27ac_MB/H3k27ac_MB_summits.bed
MB_OLIG2_summi=./Olig2_MB/Olig2_MB_summits.bed


# DeepTools:

# -R  $CB_H3K27AC_summit $CB_OLIG2_summit $MB_H3K27AC_summit $MB_OLIG2_summi \

computeMatrix scale-regions -S  $CB_OLIG2 \
                                $MB_OLIG2  \
                                $CB_H3K27AC \
                                $MB_H3K27AC \
                              -R  $CB_OLIG2_summit $MB_OLIG2_summi \
                              --beforeRegionStartLength 400 \
                              --regionBodyLength 1 \
                              --binSize 1 \
                              --afterRegionStartLength 400 \
                              --skipZeros -o matrix.mat.gz

plotHeatmap -m matrix.mat.gz -out ExampleHeatmap1.png 








