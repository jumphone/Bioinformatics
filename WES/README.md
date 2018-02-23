
##############################################
#   FASTQ to VCF             #
##############################################
Auther: Lu, Yulan 

step1: cd XXX/WES_tmp/data

step2: 建立一个新的文件夹（如17Y084）, 并将样本的R1.fastq和R2.fastq拷贝到这个新建的文件夹中（如果名字不一样，则需要改名；如果是压缩文件，则需要解压文件。即需要把名字变成R1.fastq和R2.fastq）

step3: cd /home/luyulan/project/WES_tmp

step4: sh bwa_YL.sh 新建立的文件夹名字             # demo: sh bwa_YL.sh 17Y084

       step4.1: 运行BWA                            # 把fastq中的reads匹配到参考基因组上
                                                   # -t 10 : 用10个核来跑
                                                   # 输入：R1.fastq、R2.fastq、参考基因组
                                                   # 输出：.sam文件(34G)
                                                   # 耗费40分钟
       
       step4.2: 将.sam转成.bam，并将.bam进行排序   # 输入：.sam文件
                                                   # 输出：.sort.bam文件(6G)
                                                   # 耗时38分钟
       
       step4.3: 将.sam文件去掉                     # 为了节省空间

       step4.4: 去掉PCR过程产生的duplicate         # 输入：.sort.bam文件
                                                   # 输出：.dedup.bam文件(6.1G)
                                                   # 耗时31分钟 

       step4.5: 将.sort.bam文件去掉                # 为了节省空间

       step4.6: 建立BAM文件的索引（index）         # 对.sort.bam文件建立index
                                                   # 输入：.dedup.bam文件
                                                   # 输出：.dedup.bai文件
                                                   # 耗时3分钟

step5: sh call_WES.sh 新建立的文件夹名字           # demo: sh call_WES.sh 17Y084
  
       step5.1: 对indel位点进行重新匹配（realign） # 因为indel周围出现匹配（mapping）错误的概率较高，而通过对已知的可靠indel位点（如Mills_and_1000G_gold_standard.indels.hg19.sites.vcf与1000G_phase1.indels.hg19.sites.vcf）周围进行分析，可以大大减少分析样本由于indel导致的很多假阳性snp。对indel周围进行重新匹配（realignment），需要两步操作：
         step5.1.1: 获取需要进行realign的位置信息          # 用GenomeAnalysisTK.jar中的RealignerTargetCreator
                                                           # 输入：.dedup.bam文件
                                                           # 输出：.realn.intervals文件

         step5.1.2: 对需要进行realign的位置进行realign     # 用GenomeAnalysisTK.jar中的IndelRealigner
                                                           # 输入：.realn.intervals文件、.dedup.bam文件
                                                           # 输出：.realn.bam文件
                                                           # step5.1.1与step5.1.2一共耗时73.89分钟（1.23个小时）
 
       step5.2: 对realign后得到的结果进行校正      # 因为mapping得到的分值存在偏态（不能反映真实情况），需要尽量消除偏态才能得到可靠、准确的结果。主要是通过对已知的可靠位点（如dbsnp_138.hg19.vcf、Mills_and_1000G_gold_standard.indels.hg19.sites.vcf与1000G_phase1.indels.hg19.sites.vcf）周围进行分析，对分析样本中的低覆盖度（<10x）位点进行校正，可以消除很多假阳性。
         step5.2.1: 获取需要进行校正的位置信息             # 用GenomeAnalysisTK.jar中的BaseRecalibrator
                                                           # 输入：.realn.bam文件
                                                           # 输出：.recal_date.grp文件
                                                           
         step5.2.2: 将经过校正后的数据输出到新的bam文件中  # 用GenomeAnalysisTK.jar中的PrintReads
                                                           # 输入：.recal_date.grp、.realn.bam文件
                                                           # 输出：.recal.bam文件
                                                           # 耗时3小时          
  
       step5.3: 建立BAM文件的索引（index）                 # 对.recal.bam文件建立index
                                                           # 输入：.recal.bam
                                                           # 输出：.recal.bai
                                                           # 耗时14分钟

       step5.4: 生成原始的VCF文件                          # 用GenomeAnalysisTK.jar中的HaplotypeCaller
                                                           # 输入：.recal.bam文件
                                                           # 输出：.recal.vcf文件(91M)
                                                           # 耗时9小时 
                               
       step5.5: 识别SNP                                    # 输入：.recal.vcf 文件
                                                           # 输出：.snp.vqr.vcf（106M）

       step5.6: 识别indel                                   
       step5.7: 合并snp和indel文件                         # 最终输出文件：.both.vqr.vcf文件(107M)
                                                           # step5.5、step5.6、step5.7耗时很短 
