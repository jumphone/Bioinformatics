
library(Seurat)
library(dplyr)
library(Matrix)

exp_data=read.table('magic.csv',header=T,row.names=1,sep='\t',check.names=F)
exp_data=t(exp_data)

pbmc <- CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes = 0, project = "Project")

pdf('VAR.pdf')
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()

length(x = pbmc@var.genes)

pbmc <- ScaleData(object = pbmc)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

pdf('PCA.pdf')
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
PCAPlot(object = pbmc, dim.1 = 2, dim.2 = 3)
PCAPlot(object = pbmc, dim.1 = 3, dim.2 = 4)
PCElbowPlot(object = pbmc)
dev.off()

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:5, resolution = 0.3, print.output = 0, save.SNN = TRUE,force.recalc=T)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:5, do.fast = TRUE)

pdf('TSNE.pdf')
TSNEPlot(object = pbmc,do.label = TRUE)
dev.off()


Last login: Mon Jul 16 10:07:43 on ttys000
EA18-00314:~ zha8dh$ ssh zhangfeng@202.120.224.143 -p 13579
zhangfeng@202.120.224.143's password: 
Last login: Sat Jul 14 06:50:30 2018 from 205.142.197.113
[zhangfeng@rna ~]$ ls
2016_zf_ppt.pdf                              index.html.2
CASE.zip                                     knitr.mp4
Desktop                                      listener
DhsTgene.txt.sorted.bed                      mTOR
Documents                                    mail
Downloads                                    makedir.sh
FMD_calls_commercial                         mysql
GenePANDA                                    ncbi
HK_genes.txt                                 ncbi_error_report.xml
HK_genes.txt.mouse                           ncomms1324.pdf
LAB_DATA                                     nltk_data
Lab_files                                    node1_CPU.txt
LicenseAdministration.pdf                    node1_mem.txt
NOT_FINISHED_LIST.txt                        node1_summary.txt
Nature.zip                                   nohup.out
Olig2_MB_peaks.bed                           ok.sh
PACK.zip                                     ok_joined.bed
R                                            perl5
RNA_e_data                                   phenobayes.py
Requirements for  Homework (Important).pdf   ppt
Rplots.pdf                                   project
SRR1658570_allValidPairs                     pypi-stats-cache
Sample                                       qrcode_for_gh_8341c79242cd_258.jpg
TCGA_token                                   qrcode_for_gh_8341c79242cd_430.jpg
TMP                                          question.odt
TXT2HPO.zip                                  run1663_5000_dge.txt
__MACOSX                                     run1663_5000_dge.txt.zip
abyss-1.5.2                                  run2047_lane1_read1_index1-TS=1.fa
abyss-1.5.2.tar.gz                           samtools-1.1
access_log                                   server
ana                                          sprint-0.1.7.tar.gz
bam2wig-master.zip                           sprint.zip
bed_chr_1.bed.gz                             summary_20160923.pdf
bin                                          system_bac
biostat                                      tabix-0.2.6
biostat2015                                  test
bwa-0.5.9rc1                                 testnew
dataSummary.html                             tmp
dead.letter                                  tmp.bam
disk                                         tmp.pdb
download                                     tmp.pdf
fastq_folder                                 tmp.sh
gdc-user-token.2017-07-23T01-42-52.300Z.txt  tmp____
gdc-user-token.2017-08-24T11-26-04.362Z.txt  tmp_lst.txt
gdc-user-token.2018-01-09T22_39_13.900Z.txt  tmp_macs2_model.r
gene_info                                    tmp_macs2_peaks.narrowPeak
getweb                                       tmp_macs2_peaks.xls
go                                           tmp_macs2_summits.bed
gsea-3.0.jar                                 tmptmp
gsea_home                                    tmptmptmp.txt
guide.pdf                                    tools
hgcentral.sql                                trim_galore
httpd.conf.bac                               try.pl
igv                                          workflow.odg
impress.js-master.zip                        workspace
index.html                                   www
index.html.1
[zhangfeng@rna ~]$ cd disk/project/
[zhangfeng@rna project]$ ksl
bash: ksl: command not found...
[zhangfeng@rna project]$ ls
ACTGCODE          CNV         GSE9782          HiChip      MECP2         RECT                    SSN_NA         YANZIJUN        sprint_record
ANGEL             Cuffdiff    GSEA             IGV         MEDSCIENCES   RECTK                   Single_cell    Yaqi            test
BWA               DONGa       GSVA             IR_UE       Methy         RECT_test               TCGA           cell            tmp
CAFA              DYG         GWAS             InterVar    MutationRate  REFSEQ2HGNC             TCGA_2018      ctDNA           tools
CALL_variant      ENCODE      GWAS_catlog      JINCHENG    NIPT          RNA_e_new               TCGA_2018_PUB  data            translation
CDD               GATK        GenePanda        KETI3_DATA  Natalie       RNAsplicing             TCGA_RAW       insert_GM12878  tryRECT
CHIP              GBM_MB      Genome_assembly  L1000       OCEAN         SINGLE_CELL             TCGA_old       m6A
CHIP_PEAK         GCA         HISAT2           LIGUO       OMIM          SPRINT                  TF_ENRICH      mTOR
CHPO              GESS        HPO              LUOZAILI    O_O           SPRINTX                 TTT            matlab.doc
CHUNTAO           GIREMI      HPO_DISEASE      LUQING      PersonPCC     SPRINT_7                VEP            mirRNA
CHUNTAO_20180126  GOBayes     HPO_chinese      LZX         PhenoBayes    SPRINT_bioinfo_compare  WANG           nonCDS
CMAP_MET          GO_chinese  HUXIAOHUA        MAGIC       Pubmed2HPO    SSN_CHENLUONAN          YANGNAN        old_CHIP
[zhangfeng@rna project]$ ls
ACTGCODE      CHPO              ENCODE      GSE9782          HPO          JINCHENG    MECP2         O_O          RNA_e_new               SSN_NA         TTT       data            sprint_record
ANGEL         CHUNTAO           GATK        GSEA             HPO_DISEASE  KETI3_DATA  MEDSCIENCES   PersonPCC    RNAsplicing             Single_cell    VEP       insert_GM12878  test
BWA           CHUNTAO_20180126  GBM_MB      GSVA             HPO_chinese  L1000       Methy         PhenoBayes   SINGLE_CELL             TCGA           WANG      m6A             tmp
CAFA          CMAP_MET          GCA         GWAS             HUXIAOHUA    LIGUO       MutationRate  Pubmed2HPO   SPRINT                  TCGA_2018      YANGNAN   mTOR            tools
CALL_variant  CNV               GESS        GWAS_catlog      HiChip       LUOZAILI    NIPT          RECT         SPRINTX                 TCGA_2018_PUB  YANZIJUN  matlab.doc      translation
CDD           Cuffdiff          GIREMI      GenePanda        IGV          LUQING      Natalie       RECTK        SPRINT_7                TCGA_RAW       Yaqi      mirRNA          tryRECT
CHIP          DONGa             GOBayes     Genome_assembly  IR_UE        LZX         OCEAN         RECT_test    SPRINT_bioinfo_compare  TCGA_old       cell      nonCDS
CHIP_PEAK     DYG               GO_chinese  HISAT2           InterVar     MAGIC       OMIM          REFSEQ2HGNC  SSN_CHENLUONAN          TF_ENRICH      ctDNA     old_CHIP
[zhangfeng@rna project]$ cd MAGIC/
[zhangfeng@rna MAGIC]$ ls
2Magic.py  GCA-master      GCA.zip       GCA_OUT.data.zmat             GCA_OUT.data.zmat.gca_result.saved_RData  GSE70630_OG_processed_data_v2.txt.cleaned.txt  Natalie    nohup.out  tmagic.csv
DIPG       GCA-master.zip  GCA_OUT.data  GCA_OUT.data.zmat.gca_result  GCA_OUT.data.zmat_TMP_0.534015207731      MATRIX.txt                                     magic.csv  prepare.R  transpose.R
[zhangfeng@rna MAGIC]$ cd Natalie/
[zhangfeng@rna Natalie]$ ls
1QC.R      3Seurat.R  MATRIX.txt                       PCA.pdf  SAVE.Robj    TSNE.pdf  magic.csv    nohup.out             top10.txt
2Magic.py  HEAT.pdf   MPNST4_raw_gene_bc_matrices.zip  QC.pdf   Seurat.Robj  VAR.pdf   markers.txt  raw_gene_bc_matrices
[zhangfeng@rna Natalie]$ ll
total 694804
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       417 Jul 14 04:01 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 316002197 Jul 14 04:02 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng       402 Jul 14 04:02 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ vi 2Magic.py 
[zhangfeng@rna Natalie]$ which python           
alias python='/usr/local/bin/python'
	/usr/local/bin/python
[zhangfeng@rna Natalie]$ which python3
alias python3='/home/sherry/bin/bin/python3'
	/home/sherry/bin/bin/python3
[zhangfeng@rna Natalie]$ nohup /home/sherry/bin/bin/python3 2Magic.py &
[1] 130419
[zhangfeng@rna Natalie]$ nohup: ignoring input and appending output to 'nohup.out'

[zhangfeng@rna Natalie]$ ll
total 694804
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       443 Jul 16 23:10 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 316002197 Jul 14 04:02 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng       402 Jul 14 04:02 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ ll
total 694804
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       443 Jul 16 23:10 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 316002197 Jul 14 04:02 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng       402 Jul 14 04:02 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ ll
total 694804
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       443 Jul 16 23:10 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 316002197 Jul 14 04:02 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng       402 Jul 14 04:02 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ top

top - 23:14:02 up 60 days, 10:41, 20 users,  load average: 3.92, 4.76, 4.97
Tasks: 947 total,   4 running, 940 sleeping,   1 stopped,   2 zombie
%Cpu(s):  7.5 us,  0.0 sy,  0.0 ni, 92.5 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
KiB Mem : 13161577+total, 72817120 free,  8831264 used, 49967396 buff/cache
KiB Swap:  4194300 total,  2887268 free,  1307032 used. 12218638+avail Mem 

   PID USER      PR  NI    VIRT    RES    SHR S  %CPU %MEM     TIME+ COMMAND                                                                                                                         
128335 yzj       20   0 2832656 2.319g  13480 R 100.0  1.8  58:42.15 python src/SASE_hunter.py LAML 1000                                                                                             
128378 yzj       20   0 2824720 2.312g  13480 R 100.0  1.8  58:36.95 python src/SASE_hunter.py LUAD 1000                                                                                             
128421 yzj       20   0 2826512 2.313g  13480 R  99.7  1.8  58:32.44 python src/SASE_hunter.py OV 1000                                                                                               
  2484 zhangfe+  20   0  781500  42104   2080 S   0.3  0.0 501:48.76 /home/zhangfeng/server/tianlab/tianlab_env/bin/python /home/zhangfeng/server/tianlab/manage.py runserver 9500                   
130462 zhangfe+  20   0   43504   2884   1364 R   0.3  0.0   0:00.76 top                                                                                                                             
     1 root      20   0  191532   2680   1408 S   0.0  0.0   0:32.86 /usr/lib/systemd/systemd --switched-root --system --deserialize 23                                                              
     2 root      20   0       0      0      0 S   0.0  0.0   0:01.58 [kthreadd]                                                                                                                      
     3 root      20   0       0      0      0 S   0.0  0.0   0:43.72 [ksoftirqd/0]                                                                                                                   
     5 root       0 -20       0      0      0 S   0.0  0.0   0:00.00 [kworker/0:0H]                                                                                                                  
     6 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [kworker/u384:0]                                                                                                                
     8 root      rt   0       0      0      0 S   0.0  0.0   3:22.52 [migration/0]                                                                                                                   
     9 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcu_bh]                                                                                                                        
    10 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/0]                                                                                                                       
    11 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/1]                                                                                                                       
    12 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/2]                                                                                                                       
    13 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/3]                                                                                                                       
    14 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/4]                                                                                                                       
    15 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/5]                                                                                                                       
    16 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/6]                                                                                                                       
    17 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/7]                                                                                                                       
    18 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/8]                                                                                                                       
    19 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/9]                                                                                                                       
    20 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/10]                                                                                                                      
    21 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/11]                                                                                                                      
    22 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/12]                                                                                                                      
    23 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/13]                                                                                                                      
    24 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/14]                                                                                                                      
    25 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/15]                                                                                                                      
    26 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/16]                                                                                                                      
    27 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/17]                                                                                                                      
    28 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/18]                                                                                                                      
    29 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/19]                                                                                                                      
    30 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/20]                                                                                                                      
    31 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/21]                                                                                                                      
    32 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/22]                                                                                                                      
    33 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/23]                                                                                                                      
    34 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/24]                                                                                                                      
    35 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/25]                                                                                                                      
    36 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/26]                                                                                                                      
    37 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/27]                                                                                                                      
    38 root      20   0       0      0      0 S   0.0  0.0   0:00.00 [rcuob/28]                                                                                                                      
[1]+  Done                    nohup /home/sherry/bin/bin/python3 2Magic.py
[zhangfeng@rna Natalie]$ ll
total 694684
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       443 Jul 16 23:10 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 315876052 Jul 16 23:12 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng       833 Jul 16 23:11 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ less magic.csv  
[zhangfeng@rna Natalie]$ cut -f 2 magic.csv | less
[zhangfeng@rna Natalie]$ vi 3Seurat.R 
[zhangfeng@rna Natalie]$ pwd
/home/zhangfeng/disk/project/MAGIC/Natalie
[zhangfeng@rna Natalie]$ zip -r magic.csv.zip magic.csv 
  adding: magic.csv (deflated 56%)
[zhangfeng@rna Natalie]$ ls
1QC.R      3Seurat.R  MATRIX.txt                       PCA.pdf  SAVE.Robj    TSNE.pdf  magic.csv      markers.txt  raw_gene_bc_matrices
2Magic.py  HEAT.pdf   MPNST4_raw_gene_bc_matrices.zip  QC.pdf   Seurat.Robj  VAR.pdf   magic.csv.zip  nohup.out    top10.txt
[zhangfeng@rna Natalie]$ vi 3Seurat.R 
[zhangfeng@rna Natalie]$ vi 2Magic.py 
[zhangfeng@rna Natalie]$ nohup /home/sherry/bin/bin/python3 2Magic.py &
[1] 131083
[zhangfeng@rna Natalie]$ nohup: ignoring input and appending output to 'nohup.out'

[zhangfeng@rna Natalie]$ 
[1]+  Done                    nohup /home/sherry/bin/bin/python3 2Magic.py
[zhangfeng@rna Natalie]$ 
[zhangfeng@rna Natalie]$ 
[zhangfeng@rna Natalie]$ ls
1QC.R      3Seurat.R  MATRIX.txt                       PCA.pdf  SAVE.Robj    TSNE.pdf  magic.csv      markers.txt  raw_gene_bc_matrices
2Magic.py  HEAT.pdf   MPNST4_raw_gene_bc_matrices.zip  QC.pdf   Seurat.Robj  VAR.pdf   magic.csv.zip  nohup.out    top10.txt
[zhangfeng@rna Natalie]$ less MATRIX.txt 
[zhangfeng@rna Natalie]$ zip -r MATRIX.txt.zip MATRIX
	zip warning: name not matched: MATRIX

zip error: Nothing to do! (try: zip -r MATRIX.txt.zip . -i MATRIX)
[zhangfeng@rna Natalie]$ zip -r MATRIX.txt.zip MATRIX.txt 
  adding: MATRIX.txt (deflated 96%)
[zhangfeng@rna Natalie]$ ll
total 1047148
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       439 Jul 16 23:29 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-rw-r-- 1 zhangfeng zhangfeng   1218104 Jul 16 23:31 MATRIX.txt.zip
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 315914901 Jul 16 23:29 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng 138773582 Jul 16 23:16 magic.csv.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng      1235 Jul 16 23:29 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ 
[zhangfeng@rna Natalie]$ ll
total 1047148
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       439 Jul 16 23:29 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-rw-r-- 1 zhangfeng zhangfeng   1218104 Jul 16 23:31 MATRIX.txt.zip
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 315914901 Jul 16 23:29 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng 138773582 Jul 16 23:16 magic.csv.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng      1235 Jul 16 23:29 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ top -u zhangfeng

top - 23:31:20 up 60 days, 10:58, 19 users,  load average: 3.17, 3.17, 3.70
Tasks: 941 total,   4 running, 934 sleeping,   1 stopped,   2 zombie
%Cpu(s):  7.5 us,  0.0 sy,  0.0 ni, 92.4 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
KiB Mem : 13161577+total, 72721128 free,  8788488 used, 50106164 buff/cache
KiB Swap:  4194300 total,  2888124 free,  1306176 used. 12222961+avail Mem 

   PID USER      PR  NI    VIRT    RES    SHR S  %CPU %MEM     TIME+ COMMAND                                                                                                                         
131211 zhangfe+  20   0   43372   2876   1368 R   0.8  0.0   0:00.02 top -u zhangfeng                                                                                                                
  1874 zhangfe+  20   0  270084      4      0 S   0.0  0.0   0:00.46 /home/zhangfeng/server/tianlab/tianlab_env/bin/python /home/zhangfeng/server/tianlab/manage.py runserver 9500                   
  2484 zhangfe+  20   0  781500  42104   2080 S   0.0  0.0 501:53.29 /home/zhangfeng/server/tianlab/tianlab_env/bin/python /home/zhangfeng/server/tianlab/manage.py runserver 9500                   
129998 zhangfe+  20   0  144004   2820   1524 S   0.0  0.0   0:00.18 sshd: zhangfeng@pts/17                                                                                                          
129999 zhangfe+  20   0  116700   3484   1784 S   0.0  0.0   0:00.10 -bash                                                                                                                           




































[zhangfeng@rna Natalie]$ ll
total 1047148
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       439 Jul 16 23:29 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-rw-r-- 1 zhangfeng zhangfeng   1218104 Jul 16 23:31 MATRIX.txt.zip
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 315914901 Jul 16 23:29 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng 138773582 Jul 16 23:16 magic.csv.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng      1235 Jul 16 23:29 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ less ma
ma: No such file or directory
[zhangfeng@rna Natalie]$ tail nohup.out 
Calculating PCA...
Calculated PCA in 4.86 seconds.
Calculating KNN search...
Calculated KNN search in 0.35 seconds.
Calculating affinities...
Calculated affinities in 0.65 seconds.
Calculated graph and diffusion operator in 5.97 seconds.
Calculating imputation...
Calculated imputation in 0.34 seconds.
Calculated MAGIC in 6.57 seconds.
[zhangfeng@rna Natalie]$ ll
total 1047148
-rw-rw-r-- 1 zhangfeng zhangfeng      1209 Jul 14 04:01 1QC.R
-rw-rw-r-- 1 zhangfeng zhangfeng       439 Jul 16 23:29 2Magic.py
-rw-rw-r-- 1 zhangfeng zhangfeng      2196 Jul 14 04:26 3Seurat.R
-rw-rw-r-- 1 zhangfeng zhangfeng    305894 Jul 14 04:24 HEAT.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng  30785520 Jul 14 04:01 MATRIX.txt
-rw-rw-r-- 1 zhangfeng zhangfeng   1218104 Jul 16 23:31 MATRIX.txt.zip
-rw-r--r-- 1 zhangfeng zhangfeng   9735213 Jul 14 03:53 MPNST4_raw_gene_bc_matrices.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     35089 Jul 14 04:06 PCA.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng     78847 Jul 14 04:00 QC.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 350455435 Jul 14 04:26 SAVE.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng   3793865 Jul 14 04:01 Seurat.Robj
-rw-rw-r-- 1 zhangfeng zhangfeng     18343 Jul 14 04:23 TSNE.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng    148735 Jul 14 04:06 VAR.pdf
-rw-rw-r-- 1 zhangfeng zhangfeng 315914901 Jul 16 23:29 magic.csv
-rw-rw-r-- 1 zhangfeng zhangfeng 138773582 Jul 16 23:16 magic.csv.zip
-rw-rw-r-- 1 zhangfeng zhangfeng     80102 Jul 14 04:27 markers.txt
-rw------- 1 zhangfeng zhangfeng      1235 Jul 16 23:29 nohup.out
drwxrwxr-x 3 zhangfeng zhangfeng        25 Jul 12 19:38 raw_gene_bc_matrices
-rw-rw-r-- 1 zhangfeng zhangfeng      2329 Jul 14 04:27 top10.txt
[zhangfeng@rna Natalie]$ zip -r magic.csv.zip magic.csv
updating: magic.csv^C


zip error: Interrupted (aborting)
[zhangfeng@rna Natalie]$ rm magic.csv.zip 
rm: remove regular file 'magic.csv.zip'? y
[zhangfeng@rna Natalie]$ rm MATRIX.txt.zip 
rm: remove regular file 'MATRIX.txt.zip'? y
[zhangfeng@rna Natalie]$ vi 3Seurat.R 
[zhangfeng@rna Natalie]$ vi 3Seurat.R 
[zhangfeng@rna Natalie]$ Rscript 3Seurat.R 
Loading required package: ggplot2
Loading required package: cowplot

Attaching package: 'cowplot'

The following object is masked from 'package:ggplot2':

    ggsave

Loading required package: Matrix

Attaching package: 'dplyr'

The following objects are masked from 'package:stats':

    filter, lag

The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union

^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C^C
Execution halted
[zhangfeng@rna Natalie]$ vi 3Seurat.R 
[zhangfeng@rna Natalie]$ Rscript 3Seurat.R 
Loading required package: ggplot2
Loading required package: cowplot

Attaching package: 'cowplot'

The following object is masked from 'package:ggplot2':

    ggsave

Loading required package: Matrix

Attaching package: 'dplyr'

The following objects are masked from 'package:stats':

    filter, lag

The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union

Calculating gene means
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variance to mean ratios
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Warning message:
In KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin,  :
  Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'
null device 
          1 
[1] 2114
NormalizeData has not been run, therefore ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.
ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.
Scaling data matrix
  |======================================================================| 100%
[1] "PC1"
[1] "TAF1D"    "SON"      "TRIM38"   "KCNQ1OT1" "TTC14"   
[1] ""
[1] "RPLP1"  "CNRIP1" "LGALS1" "RPL39"  "TMSB10"
[1] ""
[1] ""
[1] "PC2"
[1] "DBN1"    "IFI27L2" "SPNS1"   "SAP30BP" "STMN1"  
[1] ""
[1] "AFF1"  "MBD4"  "ZMAT3" "UCK2"  "CREM" 
[1] ""
[1] ""
[1] "PC3"
[1] "SPP1" "CTSK" "ACP5" "MMP9" "CKB" 
[1] ""
[1] "FNIP1"  "MT1E"   "NEBL"   "CTGF"   "D2HGDH"
[1] ""
[1] ""
[1] "PC4"
[1] "CADM1" "TUBG1" "DPY30" "RBM8A" "ASXL2"
[1] ""
[1] "CST4"   "ANKH"   "MT-CO3" "MXRA8"  "CD276" 
[1] ""
[1] ""
[1] "PC5"
[1] "CUEDC1" "GJA1"   "HMGB2"  "CADM1"  "PDLIM7"
[1] ""
[1] "CTSK"  "MT1E"  "CKB"   "MMP9"  "MYO1D"
[1] ""
[1] ""
null device 
          1 
null device 
          1 
# A tibble: 9 x 7
# Groups:   cluster [5]
     p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene  
     <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
1 1.04e-42    0.0413     1     1  1.07e-38 0       TMSB10
2 1.33e-36    0.0399     1     1  1.38e-32 0       LGALS1
3 6.13e-47    0.0396     1     1  6.34e-43 1       TIMP1 
4 1.45e-40    0.0442     1     1  1.50e-36 1       MT2A  
5 6.50e-26    0.0298     1     1  6.72e-22 2       MT2A  
6 1.39e-23    0.0234     1     1  1.44e-19 2       TIMP1 
7 2.62e-41    0.712      1     1  2.71e-37 3       NEAT1 
8 2.10e-38    2.03       1     1  2.18e-34 3       MALAT1
9 2.80e- 6    0.0110     1     1  2.90e- 2 4       MT2A  
null device 
          1 
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Scaling data matrix
  |======================================================================| 100%
^C
Execution halted
^C[zhangfeng@rna Natalie]$ ^C
[zhangfeng@rna Natalie]$ ^C
[zhangfeng@rna Natalie]$ Rscript 3Seurat.R 
Loading required package: ggplot2
Loading required package: cowplot

Attaching package: 'cowplot'

The following object is masked from 'package:ggplot2':

    ggsave

Loading required package: Matrix

Attaching package: 'dplyr'

The following objects are masked from 'package:stats':

    filter, lag

The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union

Calculating gene means
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variance to mean ratios
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Warning message:
In KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin,  :
  Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'
null device 
          1 
[1] 2114
NormalizeData has not been run, therefore ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.
ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.
Scaling data matrix
  |======================================================================| 100%
[1] "PC1"
[1] "TAF1D"    "SON"      "TRIM38"   "KCNQ1OT1" "TTC14"   
[1] ""
[1] "RPLP1"  "CNRIP1" "LGALS1" "RPL39"  "TMSB10"
[1] ""
[1] ""
[1] "PC2"
[1] "DBN1"    "IFI27L2" "SPNS1"   "SAP30BP" "STMN1"  
[1] ""
[1] "AFF1"  "MBD4"  "ZMAT3" "UCK2"  "CREM" 
[1] ""
[1] ""
[1] "PC3"
[1] "SPP1" "CTSK" "ACP5" "MMP9" "CKB" 
[1] ""
[1] "FNIP1"  "MT1E"   "NEBL"   "CTGF"   "D2HGDH"
[1] ""
[1] ""
[1] "PC4"
[1] "CADM1" "TUBG1" "DPY30" "RBM8A" "ASXL2"
[1] ""
[1] "CST4"   "ANKH"   "MT-CO3" "MXRA8"  "CD276" 
[1] ""
[1] ""
[1] "PC5"
[1] "CUEDC1" "GJA1"   "HMGB2"  "CADM1"  "PDLIM7"
[1] ""
[1] "CTSK"  "MT1E"  "CKB"   "MMP9"  "MYO1D"
[1] ""
[1] ""
null device 
          1 
null device 
          1 
# A tibble: 9 x 7
# Groups:   cluster [5]
     p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene  
     <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
1 1.04e-42    0.0413     1     1  1.07e-38 0       TMSB10
2 1.33e-36    0.0399     1     1  1.38e-32 0       LGALS1
3 6.13e-47    0.0396     1     1  6.34e-43 1       TIMP1 
4 1.45e-40    0.0442     1     1  1.50e-36 1       MT2A  
5 6.50e-26    0.0298     1     1  6.72e-22 2       MT2A  
6 1.39e-23    0.0234     1     1  1.44e-19 2       TIMP1 
7 2.62e-41    0.712      1     1  2.71e-37 3       NEAT1 
8 2.10e-38    2.03       1     1  2.18e-34 3       MALAT1
9 2.80e- 6    0.0110     1     1  2.90e- 2 4       MT2A  
null device 
          1 
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Scaling data matrix
  |======================================================================| 100%
# A tibble: 10 x 7
# Groups:   cluster [5]
      p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene  
      <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
 1 6.32e-10     0.593 0.354 0.187  6.54e- 6 0       SPP1  
 2 2.79e- 6     0.368 0.307 0.186  2.89e- 2 0       CTSK  
 3 3.35e- 6     0.297 0.91  0.867  3.47e- 2 1       FN1   
 4 1.22e- 4     0.308 0.91  0.88   1.00e+ 0 1       TIMP1 
 5 3.36e- 5     0.160 0.992 0.966  3.48e- 1 2       MT2A  
 6 9.52e- 3     0.113 0.407 0.581  1.00e+ 0 2       RPS29 
 7 1.02e-76     1.78  1     0.911  1.05e-72 3       MALAT1
 8 6.78e-65     1.69  0.951 0.596  7.01e-61 3       NEAT1 
 9 7.50e-23     0.345 0.876 0.611  7.76e-19 4       NEAT1 
10 1.42e- 3     0.380 0.323 0.223  1.00e+ 0 4       PSMC5 
null device 
          1 
[zhangfeng@rna Natalie]$ vi 3Seurat.R 
[zhangfeng@rna Natalie]$ vi 3Seurat.R 


pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:5, resolution = 0.3, print.output = 0, save.SNN = TRUE,force.recalc=T)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:5, do.fast = TRUE)

pdf('TSNE.pdf')
TSNEPlot(object = pbmc,do.label = TRUE)
dev.off()


#save(pbmc,file='SAVE.Robj')

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,test.use='t', logfc.threshold=0.01)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf('HEAT.pdf',width=15,height=15)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
write.table(top10,file='top10.txt',sep='\t',quote=F,row.names=T,col.names=T)
write.table(pbmc.markers,file='markers.txt',sep='\t',quote=F,row.names=T,col.names=T)
#####################


ori_exp_data=read.table('MATRIX.txt',header=T,row.names=1,sep=',',check.names=F)
ori_exp_data=t(ori_exp_data)

ori_pbmc <- CreateSeuratObject(raw.data = ori_exp_data, min.cells = 0, min.genes = 0, project = "Origin")
ori_pbmc <- NormalizeData(object = ori_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
ori_pbmc <- ScaleData(object = ori_pbmc)
ori_pbmc@ident=pbmc@ident

pbmc.markers <- FindAllMarkers(object = ori_pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,test.use='t', logfc.threshold=0.01)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
pdf('ori_HEAT.pdf',width=15,height=15)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
write.table(top10,file='ori_top10.txt',sep='\t',quote=F,row.names=T,col.names=T)
write.table(pbmc.markers,file='ori_markers.txt',sep='\t',quote=F,row.names=T,col.names=T)


###########################



