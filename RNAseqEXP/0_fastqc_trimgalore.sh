read1=A1.fastq
read2=A2.fastq
OUTPUT=./

FASTQC=/home/zhangfeng/tools/FastQC/fastqc/fastqc
TRIMGALORE=/home/zhangfeng/tools/TrimGalore-0.4.4/trim_galore

$FASTQC $read1 $read2
$TRIMGALORE --illumina --paired  $read1 $read2 -o $OUTPUT
