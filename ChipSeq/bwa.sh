cpu=$1
refgenome=$2
treat=$3
bwa='/home/zhangfeng/tools/bwa-0.7.12/bwa'
samtools='/home/zhangfeng/disk/project/SPRINT_bioinfo_compare/softwares/samtools-1.4/samtools'
$bwa mem  -M -t $cpu $refgenome $treat > $treat\.sam
$samtools view -bS $treat\.sam > $treat\.bam
$samtools sort $treat\.bam > $treat\.sorted.bam
$samtools rmdup $treat\.sorted.bam $treat\.sorted.dedup.bam
