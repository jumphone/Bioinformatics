CPU=10
BISMARK='XXX/Methy/tools/Bismark_v0.19.0/bismark --score_min L,0,-0.2 -I 0 -X 1000 '
REFGENOME=XXX/Methy/reference/hg19
OUT=XXX/Methy/data/SAM
DATA=XXX/Methy/data/
TRIM=XXX/tools/trim_galore


TAG=AO2
READ1=$DATA$TAG\_1.fastq.gz
READ2=$DATA$TAG\_2.fastq.gz
$TRIM --paired $READ1 $READ2 -o $DATA
$BISMARK -p $CPU $REFGENOME -1 $DATA$TAG\_1_val_1.fq.gz -2 $DATA$TAG\_2_val_2.fq.gz --sam -o $OUT
