B=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/run1978_lane12_Pdgfra-GFP/FILTERED
M=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/REF/mm10_rmsk.gtf
GTF=/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf
BAM=/home/disk/RNAediting_Cancer/tmp_single_lu_drop/run1978_lane12_Pdgfra-GFP/run1978_lane12_read1_indexN709=Drop-seq_Redo_Potter_Lu_Pdgfrb_gfp_2_13_17_N709.fastq_BCnum_3184/picard.bam.merged.bam

PATH=''
export PATH='/home/disk/RNAediting_Cancer/tmp_single_lu_drop/run1978_lane12_Pdgfra-GFP/run1978_lane12_read1_indexN709=Drop-seq_Redo_Potter_Lu_Pdgfrb_gfp_2_13_17_N709.fastq_BCnum_3184/download/samtools-1.7/'
export PATH=$PATH':/usr/bin/'

/usr/local/bin/velocyto run  -b $B -o $BAM\_velocyto -m $M  $BAM $GTF
