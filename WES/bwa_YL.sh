tag=$1
dir=/home/luyulan/project/WES_tmp/data
dir_out=/home/luyulan/project/WES_tmp/data/${tag}

picard_dir=/home/luyulan/Tools/picard-tools-1.131
gatk_dir=/home/luyulan/Tools/GATK/GenomeAnalysisTK-3.3-0
known_dir=${gatk_dir}/bundle
ref_genome=${known_dir}/ucsc.hg19.fasta

fastq_p1=${dir}/${tag}/R1.fastq
fastq_p2=${dir}/${tag}/R2.fastq
input_bam=${dir_out}/${tag}

mkdir $dir_out

run_log=${dir_out}/${tag}.log

rm -rf $run_log
echo Job start on  >>$run_log
date >>$run_log
echo Input Fastq $fastq_p1 >>$run_log
echo Ref_Genome $ref_genome >>$run_log
echo Out_put $input_bam >>$run_log

#	step 0: index ref_genome

#bwa index -a bwtsw $ref_genome

#	step 1: RAW mapping
#	part a: mapping reads to .sam file

<<!EOF
echo BWA mapping start on  >>$run_log
date >>$run_log
bwa mem -t 10 -M -R '@RG\tID:'$tag'\tSM:'$tag'\tPL:illumina\tLB:lib1\tPU:unit1' $ref_genome $fastq_p1 $fastq_p2 >${input_bam}.sam
!EOF

#java -Xmx4g -jar ${picard_dir}/AddOrReplaceReadGroups.jar INPUT=${input_bam}.bam OUTPUT=${input_bam}.RG.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${tag}

<<!EOF
#	part b:	transfer .sam to .bam and sort
echo Transfer .sam to .bam on  >>$run_log
date >>$run_log
java -Xmx4g -jar ${picard_dir}/SortSam.jar INPUT= ${input_bam}.sam OUTPUT= ${input_bam}.sort.bam SORT_ORDER= coordinate
!EOF

<<!EOF
wait
rm -rf ${input_bam}.sam 
!EOF

<<!EOF
# part c: filter PCR duplicates
echo Filter PCR dup on  >>$run_log
date >>$run_log
java -Xmx4g -jar ${picard_dir}/MarkDuplicates.jar INPUT= ${input_bam}.sort.bam OUTPUT= ${input_bam}.dedup.bam METRICS_FILE= ${input_bam}.dedup.metrics 
!EOF

<<!EOF
wait
rm -rf ${input_bam}.sort.bam
!EOF

#<<!EOF
java -Xmx4g -jar ${picard_dir}/BuildBamIndex.jar INPUT= ${input_bam}.dedup.bam
#!EOF 
