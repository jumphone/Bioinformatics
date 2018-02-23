cpu=$1
refgenome=$2
treat=$3
bwa='/home/zhangfeng/tools/bwa-0.7.12/bwa'
samtools='/home/zhangfeng/tools/samtools-1.2/samtools'
$bwa mem  -M -t $cpu $refgenome $treat > $treat\.sam
$samtools view -bS $treat\.sam > $treat\.bam
$samtools sort $treat\.bam $treat\.bam.sorted
$samtools rmdup $treat\.bam.sorted.bam $treat\.bam.sorted.dedup.bam
