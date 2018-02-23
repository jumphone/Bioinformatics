#WP=/Users/jumphone/Desktop/CCHMC_LAB/GSEA/Yaqi_20180119/
EXP=$1 #$WP\DIPG_IV_vs_NSC.txt
CLS=$2 #$WP\DIPG_IV_vs_NSC.cls
GMT=$3 #$WP\GENESET/c2.cp.kegg.v6.1.symbols_KEGG.gmt
OUT=$4
NUM=2000
SET_MIN=5
SET_MAX=1000
SEED=13579


java -Xmx2048m -cp /Users/jumphone/Desktop/CCHMC_LAB/Projects/GSEA/gsea3/ xtools.gsea.Gsea \
-out $OUT \
-rnd_seed $SEED \
-set_min $SET_MIN \
-set_max $SET_MAX \
-res $EXP \
-cls $CLS \
-gmx $GMT \
-plot_top_x $NUM \
-permute gene_set \
-collapse false \
