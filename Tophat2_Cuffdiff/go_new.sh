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

mm10_gft='/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf'
mm9_gft='/home/genomewide/annotation/mm9/Mus_musculus.NCBIM37.67.chr.gtf'
mm10_index='/home/genomewide/refgenome/mm10/mm10'
mm10_fa='/home/genomewide/refgenome/mm10/mm10.fa'

wp=/home/disk/lynn/data/Cindy_RNA_seq/FASTQ_03062018/

read1=$wp\1052_S71_L007_R1_001.fastq
read2=$wp\1052_S71_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=$wp\1054_S72_L007_R1_001.fastq
read2=$wp\1054_S72_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=$wp\1055_S73_L007_R1_001.fastq
read2=$wp\1055_S73_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=$wp\1061_S74_L007_R1_001.fastq
read2=$wp\1061_S74_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=$wp\1078_S75_L007_R1_001.fastq
read2=$wp\1078_S75_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=$wp\1081_S76_L007_R1_001.fastq
read2=$wp\1081_S76_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=$wp\1089_S77_L007_R1_001.fastq
read2=$wp\1089_S77_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=$wp\1090_S78_L007_R1_001.fastq
read2=$wp\1090_S78_L007_R2_001.fastq
tophat2 -p 8 -o $read1\.tophatdir $mm10_index $read1 $read2
$cufflinks -p 8 -u -G $mm10_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

GTF=/home/disk/lynn/data/Cindy_RNA_seq/FASTQ_03062018/COMBINED_GTF
ASM=/home/disk/lynn/data/Cindy_RNA_seq/FASTQ_03062018/ASM.txt

$cuffmerge -g $mm10_gtf -s $mm10_fa -p 10  -o $GTF  $ASM

$KRAS_1052=$wp\1052_S71_L007_R1_001.fastq.tophatdir/accepted_hits.bam
$KRAS_1054=$wp\1054_S72_L007_R1_001.fastq.tophatdir/accepted_hits.bam
$KRAS_1078=$wp\1078_S75_L007_R1_001.fastq.tophatdir/accepted_hits.bam
$KRAS_1089=$wp\1089_S77_L007_R1_001.fastq.tophatdir/accepted_hits.bam
$CONTROL_1055=$wp\1055_S73_L007_R1_001.fastq.tophatdir/accepted_hits.bam
$CONTROL_1061=$wp\1061_S74_L007_R1_001.fastq.tophatdir/accepted_hits.bam
$CONTROL_1081=$wp\1081_S76_L007_R1_001.fastq.tophatdir/accepted_hits.bam
$CONTROL_1090=$wp\1090_S78_L007_R1_001.fastq.tophatdir/accepted_hits.bam

DIF=/home/disk/lynn/data/Cindy_RNA_seq/FASTQ_03062018/DIF
$cuffdiff -o $DIF  -p 10 -L KRAS,CONTROL -u $GTF\/merged.gtf $KRAS_1052,$KRAS_1054,$KRAS_1078,$KRAS_1089 $CONTROL_1055,$CONTROL_1061,$CONTROL_1090



