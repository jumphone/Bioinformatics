#https://github.com/nservant/HiC-Pro


#/home/disk/zhangfeng/project/HiChip/download/HiC-Pro-bin/HiC-Pro_2.11.1/bin/utils/digest_genome.py -r mboi -o MBOI_resfrag_rn5.bed rn5.fa

HicPro=/home/disk/zhangfeng/project/HiChip/download/HiC-Pro-bin/HiC-Pro_2.11.1/bin/HiC-Pro
INPUT=/home/disk/zhangfeng/project/HiChip/rawdata
OUTPUT=/home/disk/zhangfeng/project/HiChip/output
CONFIG=/home/disk/zhangfeng/project/HiChip/CONFIG1.txt

$HicPro -i $INPUT -o $OUTPUT -c $CONFIG

