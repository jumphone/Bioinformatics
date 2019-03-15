


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







