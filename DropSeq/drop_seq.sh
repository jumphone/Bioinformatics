read1=$1
read2=$2
BCnum=$3
CPU=10
tmp=$read1\_BCnum_$BCnum\/
bam=$tmp\picard.bam

refgenome_seq=/home/genomewide/refgenome/mm10/mm10.fa
refgenome_star=/home/genomewide/refgenome/mm10/mm10_star
annotation_gtf=/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf
annotation_stat=/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf.bed.stat
picard_tmp_dir=/home/RNAediting/picardtools_tmp

#Tools
picard=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar
tagbam=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended
filterbam=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/FilterBAM
trimstarting=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/TrimStartingSequence
polyATrimmer=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/PolyATrimmer
star=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/STAR-2.5.3a/bin/Linux_x86_64/STAR
TagReadWithGeneExon=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/TagReadWithGeneExon
DetectBeadSynthesisErrors=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/DetectBeadSynthesisErrors
DigitalExpression=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/tools/Drop-seq_tools-1.12/DigitalExpression

#Drop-seq pipeline

#Make tmp dir
mkdir $tmp
#fq2bam
java -jar $picard FastqToSam F1=$read1 F2=$read2 O=$tmp\picard.bam SO=queryname SAMPLE_NAME=drop_seq TMP_DIR=$picard_tmp_dir
#Add_tag
$tagbam I=$bam O=$bam\.tagged.cell.bam BARCODED_READ=1 BASE_RANGE=1-12 TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 BASE_QUALITY=10 SUMMARY=$bam\.tag_summary
$tagbam I=$bam\.tagged.cell.bam O=$bam\.tagged.cellmol.bam BARCODED_READ=1 BASE_RANGE=13-20 TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 BASE_QUALITY=10 SUMMARY=$bam\.tag_summary_cellmol DISCARD_READ=True
#Filter_bam
$filterbam TAG_REJECT=XQ INPUT=$bam\.tagged.cellmol.bam OUTPUT=$bam\.filtered.bam
#QC
$trimstarting INPUT=$bam\.filtered.bam OUTPUT=$bam\.filtered_trimmed.bam SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5
$polyATrimmer INPUT=$bam\.filtered_trimmed.bam OUTPUT=$bam\.filtered_trimmed_polyA.bam MISMATCHES=0 NUM_BASES=6
#Alignment
java -Xmx4g -jar $picard SamToFastq INPUT=$bam\.filtered_trimmed_polyA.bam FASTQ=$bam\.filtered_trimmed_polyA.fastq TMP_DIR=$picard_tmp_dir
$star --genomeDir $refgenome_star --readFilesIn $bam\.filtered_trimmed_polyA.fastq --outFileNamePrefix $bam\.filtered_trimmed_polyA.fastq.star --runThreadN $CPU
java -Xmx4g -jar $picard SortSam I=$bam\.filtered_trimmed_polyA.fastq.starAligned.out.sam O=$bam\.filtered_trimmed_polyA.fastq.sorted.bam SO=queryname  TMP_DIR=$picard_tmp_dir
java -Xmx4g -jar $picard MergeBamAlignment REFERENCE_SEQUENCE=$refgenome_seq UNMAPPED_BAM=$bam\.filtered_trimmed_polyA.bam ALIGNED_BAM=$bam\.filtered_trimmed_polyA.fastq.sorted.bam  OUTPUT=$bam\.merged.bam  PAIRED_RUN=false VALIDATION_STRINGENCY=LENIENT INCLUDE_SECONDARY_ALIGNMENTS=false TMP_DIR=$picard_tmp_dir
#Quantification
$TagReadWithGeneExon I=$bam\.merged.bam O=$bam\.merged_exon_tagged.bam ANNOTATIONS_FILE=$annotation_gtf TAG=GE TMP_DIR=$picard_tmp_dir
$DetectBeadSynthesisErrors I=$bam\.merged_exon_tagged.bam  O=$bam\.clean.bam OUTPUT_STATS=$bam\.clean.bam.stat SUMMARY=$bam\.clean.bam.summary NUM_BARCODES=$BCnum PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC TMP_DIR=$picard_tmp_dir
$DigitalExpression I=$bam\.clean.bam O=$bam\.clean.bam.dge.txt SUMMARY=$bam\.clean.bam.dge.summary NUM_CORE_BARCODES=$BCnum TMP_DIR=$picard_tmp_dir

