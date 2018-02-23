tag=$1
dir=/home/luyulan/project/WES_tmp/data

picard_dir=/home/luyulan/Tools/picard-tools-1.131
gatk_dir=/home/luyulan/Tools/GATK/GenomeAnalysisTK-3.3-0

dir_out=/home/luyulan/project/WES_tmp/data/${tag}

known_dir=${gatk_dir}/bundle
ref_genome=${known_dir}/ucsc.hg19.fasta

fastq_p1=${dir}/${tag}/R1.fastq
fastq_p2=${dir}/${tag}/R2.fastq
input_bam=${dir_out}/${tag}
mkdir ${dir}/${tag}
run_log=${dir}/${tag}/${tag}.log


rm -rf $run_log
echo Job start on  >>$run_log
date >>$run_log
echo Ref_Genome $ref_genome >>$run_log
echo Out_put $input_bam >>$run_log

#	step 2: local realignment for INDEL
<<!EOF
#	part a: define realignment region with known SNP (in vcf)

echo Define realignment region with known SNP on  >>$run_log
date >>$run_log
java -Xmx4g -jar ${gatk_dir}/GenomeAnalysisTK.jar -R $ref_genome -nt 10 -T RealignerTargetCreator	-I ${input_bam}.dedup.bam	-o ${input_bam}.realn.intervals -known ${known_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known ${known_dir}/1000G_phase1.indels.hg19.sites.vcf 
#-fixMisencodedQuals 
#--fix_misencoded_quality_scores -fixMisencodedQuals

# part b: realign in selected region

echo Realign in selected region on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar -R $ref_genome -T IndelRealigner -targetIntervals ${input_bam}.realn.intervals -I ${input_bam}.dedup.bam -o ${input_bam}.realn.bam -known ${known_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known ${known_dir}/1000G_phase1.indels.hg19.sites.vcf
#-fixMisencodedQuals 
#--fix_misencoded_quality_scores -fixMisencodedQuals
!EOF

#        step 3: base quality score recalibration
<<!EOF
#        part a: recalibration-file based on known data

echo Recalibrate file on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar -R $ref_genome -nct 10 -T BaseRecalibrator -I ${input_bam}.realn.bam -knownSites ${known_dir}/dbsnp_138.hg19.vcf -knownSites ${known_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites ${known_dir}/1000G_phase1.indels.hg19.sites.vcf -o ${input_bam}.recal_date.grp; #only be used fo GATK, not GATK Lite

#        part b: recalibration with previous .grp file

echo Recalibration with GRP on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar -R $ref_genome -T PrintReads -I ${input_bam}.realn.bam -BQSR ${input_bam}.recal_date.grp -o ${input_bam}.recal.bam
!EOF

#         step 4: variant calling
<<!EOF
#	  by HaplotypeCaller

echo Index recal.bam  >>$run_log
date >>$run_log

java -Xmx4g -jar ${picard_dir}/BuildBamIndex.jar INPUT= ${input_bam}.recal.bam
!EOF

<<!EOF
echo Call vcf by HaplotypeCaller on  >>$run_log
date >>$run_log

java -jar ${gatk_dir}/GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R $ref_genome -I ${input_bam}.recal.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o ${dir_out}.recal.vcf ;	# only for full version
!EOF

#	step 5: hard filter and Variant Quality Recalibration

#	part a: compare between multiple --TStranche cutoffs for SNP calling

echo Hard filter and VAQR on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar -R $ref_genome --maxGaussians 4 -T VariantRecalibrator -mode SNP -input ${dir_out}.recal.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${known_dir}/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 ${known_dir}/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${known_dir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${known_dir}/dbsnp_138.hg19.vcf -an DP -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an QD -recalFile ${dir_out}.snp.vqr.recal -tranchesFile ${dir_out}.snp.vqr.tranches -rscriptFile ${dir_out}.snp.vqr.R -nt 4 -tranche 95.0 -tranche 97.0	-tranche 99.0	-tranche 99.9	-tranche 100.0

# part b: re-VAQR based on the best --TStranche for SNP calling

echo re-VAQR on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar -R $ref_genome	--maxGaussians 2 -T VariantRecalibrator	-mode SNP	-input ${dir_out}.recal.vcf	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${known_dir}/hapmap_3.3.hg19.sites.vcf	-resource:omni,known=false,training=true,truth=true,prior=12.0 ${known_dir}/1000G_omni2.5.hg19.sites.vcf	-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${known_dir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${known_dir}/dbsnp_138.hg19.vcf	-an MQRankSum	--TStranche 99.0	-recalFile ${dir_out}.snp.vqr.recal	-tranchesFile ${dir_out}.snp.vqr.tranches	-rscriptFile ${dir_out}.snp.vqr.R	-nt 4 --minNumBadVariants 10000

# part c: apply VQSR to SNP calling

echo Apply VQSR on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar -R $ref_genome	-T ApplyRecalibration	-mode SNP	-input ${dir_out}.recal.vcf	-tranchesFile ${dir_out}.snp.vqr.tranches	-recalFile ${dir_out}.snp.vqr.recal	-o ${dir_out}.snp.vqr.vcf	--ts_filter_level 99.0

# part d: VQSR for INDEL

echo VQSR for INDEL on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar -T VariantRecalibrator	-R $ref_genome --minNumBadVariants 20000 -mode INDEL --maxGaussians 2	-std 10.0	-input ${dir_out}.snp.vqr.vcf	-resource:mills,known=true,training=true,truth=true,prior=12.0 ${known_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -an DP -an MQRankSum	-an ReadPosRankSum -an FS	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 97.0 -tranche 95.0	-recalFile ${dir_out}.indel.vqr.recal	-tranchesFile ${dir_out}.indel.vqr.tranche	-rscriptFile ${dir_out}.indel.vqr.R	

echo Apply VQSR for INDEL on  >>$run_log
date >>$run_log
java -jar ${gatk_dir}/GenomeAnalysisTK.jar	-T ApplyRecalibration	-R $ref_genome	-mode INDEL	-input ${dir_out}.snp.vqr.vcf	--ts_filter_level 99.0	-recalFile ${dir_out}.indel.vqr.recal	-tranchesFile ${dir_out}.indel.vqr.tranche	-o ${dir_out}.both.vqr.vcf
#!EOF
