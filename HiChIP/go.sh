ValidPairs=HiC_data/Human/out/hic_results/data/olig2-DIPG/olig2-DIPG_allValidPairs
GSIZE=/home/zhangfeng/disk/project/HiChip/hg19.chrom.sizes
JUICEBOXJAR=juicer_tools.1.8.9_jcuda.0.8.jar

./HiC-Pro-master/bin/utils/hicpro2juicebox.sh -i $ValidPairs  -g $GSIZE -j $JUICEBOXJAR -o OUT_TMP
INPUT=OUT_TMP/olig2-DIPG_allValidPairs.hic
java -jar juicer_tools.1.8.9_jcuda.0.8.jar hiccups  --ignore_sparsity -m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 $INPUT $INPUT\.HiCCUPS_output
