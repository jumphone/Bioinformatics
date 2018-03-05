export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64
export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/bowtie2
export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/tophat-2.1.0.Linux_x86_64
export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/samtools-1.2

cuffdiff='/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff'
cuffmerge='/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64/cuffmerge'
cufflinks='/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks'
pc_gtf='/home/zhangfeng/disk/project/data/annotation/mm9.gtf'
downloaddir='/home/zhangfeng/disk/project/data/download/'

human_gtf='/home/genomewide/annotation/hg19/Homo_sapiens.GRCh37.75.chr.gtf'



read1=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AC_1.fastq
read2=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AC_2.fastq
tophat2 -p 12 -o $read1\.tophatdir ../index/hg19 $read1 $read2
$cufflinks -p 8 -u -G $human_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AS_1.fastq
read2=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AS_2.fastq
tophat2 -p 8 -o $read1\.tophatdir ../index/hg19 $read1 $read2
$cufflinks -p 8 -u -G $human_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

$cuffmerge -g $human_gtf -s ../index/hg19.fa -p 10  -o /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/GTF  /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/assemblies.txt
AC='/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AC_1.fastq.tophatdir/accepted_hits.bam'
AS='/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AS_1.fastq.tophatdir/accepted_hits.bam'
$cuffdiff -o /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/Diff -p 10 -L AC,AS -u /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/GTF/merged.gtf $AC $AS
