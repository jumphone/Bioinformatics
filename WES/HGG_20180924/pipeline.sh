module load trimgalore/0.4.2
module load python/2.7.5
module load bwa/0.7.12
module load gatk/3.7
module load samtools/1.3
module load bcftools/1.3

TAG=$1 #HGG1-normal-brain
tag=$2 #HGG1normalbrain


OUTPUT=/users/zha8dh/tianlab/HGG/data/$TAG\/
mkdir $OUTPUT
HG19=/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
HG19_FAI=/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

read1='/users/zha8dh/tianlab/HGG/HGG_pairs_project_Nada/FASTQ/'$TAG'_R1.fastq.gz'
read2='/users/zha8dh/tianlab/HGG/HGG_pairs_project_Nada/FASTQ/'$TAG'_R2.fastq.gz'
trim_galore  --illumina --paired  $read1 $read2 -o $OUTPUT

read1=$OUTPUT''$TAG'_R1_val_1.fq.gz'
read2=$OUTPUT''$TAG'_R2_val_2.fq.gz'
bwa mem -t 10 -M -R '@RG\tID:'$tag'\tSM:'$tag'\tPL:illumina\tLB:lib1\tPU:unit1' $HG19  $read1 $read2 > $OUTPUT$TAG\.sam

PICARD=/usr/local/picard/1.89/jar/
java -Xmx4g -jar $PICARD\SortSam.jar INPUT= $OUTPUT$TAG\.sam OUTPUT= $OUTPUT$TAG\.sort.bam SORT_ORDER= coordinate
java -Xmx4g -jar $PICARD\/MarkDuplicates.jar INPUT= $OUTPUT$TAG\.sort.bam OUTPUT= $OUTPUT$TAG\.dedup.bam METRICS_FILE= $OUTPUT$TAG\.dedup.metrics
java -Xmx4g -jar $PICARD\/BuildBamIndex.jar INPUT= $OUTPUT$TAG\.dedup.bam

GATK=/usr/local/GATK/3.7/bin/
ANNO=/data/tianlab/zhangfeng/annotation/GATK/
java -Xmx4g -jar $GATK\GenomeAnalysisTK.jar -R $HG19_FAI -nt 10 -T RealignerTargetCreator       -I $OUTPUT$TAG\.dedup.bam       -o $OUTPUT$TAG\.realn.intervals -known $ANNO\/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known $ANNO\/1000G_phase1.indels.hg19.sites.vcf

java -jar $GATK\/GenomeAnalysisTK.jar -R $HG19_FAI -T IndelRealigner -targetIntervals $OUTPUT$TAG\.realn.intervals -I $OUTPUT$TAG\.dedup.bam -o $OUTPUT$TAG.realn.bam -known $ANNO\/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known $ANNO\/1000G_phase1.indels.hg19.sites.vcf

java -jar $GATK\/GenomeAnalysisTK.jar -R $HG19_FAI -nct 10 -T BaseRecalibrator -I $OUTPUT$TAG\.realn.bam -knownSites $ANNO\/dbsnp_138.hg19.vcf -knownSites $ANNO\/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $ANNO\/1000G_phase1.indels.hg19.sites.vcf -o $OUTPUT$TAG\.recal_date.grp; #only be used fo GATK, not GATK Lite

java -jar $GATK\/GenomeAnalysisTK.jar -R $HG19_FAI -T PrintReads -I $OUTPUT$TAG\.realn.bam -BQSR $OUTPUT$TAG\.recal_date.grp -o $OUTPUT$TAG\.recal.bam

java -Xmx4g -jar $PICARD\/BuildBamIndex.jar INPUT= $OUTPUT$TAG\.recal.bam

java -jar $GATK\/GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R $HG19_FAI -I $OUTPUT$TAG\.recal.bam --genotyping_mode DISCOVERY  -stand_call_conf 30 -o $OUTPUT$TAG\.recal.vcf ;     # only for full version


rm -rf $OUTPUT$TAG\.sam
rm -rf $OUTPUT$TAG\.sort.bam
rm -rf $OUTPUT$TAG\.realn.bam


java -jar $GATK\/GenomeAnalysisTK.jar -R $HG19_FAI --maxGaussians 4 -T VariantRecalibrator -mode SNP -input $OUTPUT$TAG\.recal.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ANNO\/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 $ANNO\/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 $ANNO\/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ANNO\/dbsnp_138.hg19.vcf -an DP -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an QD -recalFile $OUTPUT$TAG\.snp.vqr.recal -tranchesFile $OUTPUT$TAG\.snp.vqr.tranches -rscriptFile $OUTPUT$TAG\.snp.vqr.R -nt 4 -tranche 95.0 -tranche 97.0 -tranche 99.0   -tranche 99.9   -tranche 100.0

java -jar $GATK\/GenomeAnalysisTK.jar -R $HG19_FAI      -T ApplyRecalibration   -mode SNP       -input $OUTPUT$TAG\.recal.vcf   -tranchesFile $OUTPUT$TAG\.snp.vqr.tranches     -recalFile $OUTPUT$TAG\.snp.vqr.recal   -o $OUTPUT$TAG\.snp.vqr.vcf     --ts_filter_level 99.0

java -jar $GATK\/GenomeAnalysisTK.jar -T VariantRecalibrator    -R $HG19_FAI --minNumBadVariants 20000 -mode INDEL --maxGaussians 2     -std 10.0       -input $OUTPUT$TAG\.snp.vqr.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 $ANNO\/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -an DP -an MQRankSum  -an ReadPosRankSum -an FS       -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 97.0 -tranche 95.0  -recalFile $OUTPUT$TAG\.indel.vqr.recal -tranchesFile $OUTPUT$TAG\.indel.vqr.tranche    -rscriptFile $OUTPUT$TAG\.indel.vqr.R

java -jar $GATK\/GenomeAnalysisTK.jar   -T ApplyRecalibration   -R $HG19_FAI    -mode INDEL     -input $OUTPUT$TAG\.snp.vqr.vcf --ts_filter_level 99.0  -recalFile $OUTPUT$TAG\.indel.vqr.recal -tranchesFile $OUTPUT$TAG\.indel.vqr.tranche    -o $OUTPUT$TAG\.both.vqr.vcf

cp $OUTPUT$TAG\.both.vqr.vcf $OUTPUT\Final.gatk.vcf



OUTPUT=/users/zha8dh/tianlab/HGG/data/$TAG\/
HG19=/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
HG19_FAI=/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa


samtools mpileup -B -C50 -q 0 -t DP,SP,AD -d 8000 -ugf $HG19_FAI $OUTPUT$TAG\.recal.bam | bcftools view -v snps,indels  - > $OUTPUT\Final.samtools_all.vcf
freebayes -f $HG19_FAI $OUTPUT$TAG\.recal.bam > $OUTPUT\Final.freebayes.vcf

python2.7 /users/zha8dh/tianlab/HGG/pureVCF.py $OUTPUT\Final.gatk.vcf
python2.7 /users/zha8dh/tianlab/HGG/pureVCF.py $OUTPUT\Final.samtools_all.vcf
python2.7 /users/zha8dh/tianlab/HGG/pureVCF.py $OUTPUT\Final.freebayes.vcf

python2.7 /users/zha8dh/tianlab/HGG/combineVCF.py $OUTPUT\Final.samtools_all.vcf.pure $OUTPUT\Final.freebayes.vcf.pure $OUTPUT\Final.gatk.vcf.pure $OUTPUT\Final.combine.vcf



convert2annovar.pl -format vcf4 $OUTPUT\Final.combine.vcf > $OUTPUT\Final.combine.vcf.avinput
table_annovar.pl $OUTPUT\Final.combine.vcf.avinput /users/zha8dh/tianlab/annotation/ANNOVAR/  -buildver hg19  -protocol refGene,exac03,1000g2014oct_all,cosmic68 -operation g,f,f,f -nastring .







