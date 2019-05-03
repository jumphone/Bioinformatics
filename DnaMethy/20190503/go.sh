INPUT=DIPG4_1_bismark_bt2_pe.bedGraph
bedtools window -w 2000 -a Homo_sapiens.GRCh38.87.chr.gtf.combined.pc.bed -b $INPUT > $INPUT\.gene

INPUT=DIPG13_1_bismark_bt2_pe.bedGraph
bedtools window -w 2000 -a Homo_sapiens.GRCh38.87.chr.gtf.combined.pc.bed -b $INPUT > $INPUT\.gene

INPUT=DIPGC1_1_bismark_bt2_pe.bedGraph
bedtools window -w 2000 -a Homo_sapiens.GRCh38.87.chr.gtf.combined.pc.bed -b $INPUT > $INPUT\.gene

INPUT=DIPGC2_1_bismark_bt2_pe.bedGraph
bedtools window -w 2000 -a Homo_sapiens.GRCh38.87.chr.gtf.combined.pc.bed -b $INPUT > $INPUT\.gene

python getGene.py DIPG4_1_bismark_bt2_pe.bedGraph.gene
python getGene.py DIPG13_1_bismark_bt2_pe.bedGraph.gene
python getGene.py DIPGC1_1_bismark_bt2_pe.bedGraph.gene
python getGene.py DIPGC2_1_bismark_bt2_pe.bedGraph.gene 

Rscript combine.R