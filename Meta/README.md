转自： https://blog.csdn.net/zhouxin518/article/details/80392351

1. Illumina Hiseq PE150/250 测序

2. fastx 

3. Seqprep (https://github.com/jstjohn/SeqPrep) 和Sickle (https://github.com/najoshi/sickle) 进行质控后数据统计，基于原始测序数据，使用相应软件对其进行数据质控，剪切掉数据中的低质量及含N的reads，获得后续分析需要的高质量序列。

4. BWA去宿主后数据统计，去除宿主污染: Plants, Solanum_lycopersicum

5. Multiple_Megahit (https://github.com/voutcn/megahit) 最短contig长度 ≥ 300 bp 拼接组装与基因预测。通过相应的拼接软件，选择拼接效果最佳的序列，对结果进行ORF预测。选择核酸长度大于等于100bp的基因，并将其翻译为氨基酸序列。

6. CD-HIT 基因序列聚类相似度（Identity）≥ 0.95 基因序列聚类覆盖度 （Coverage）≥ 0.9。通过CD-HIT软件对样本预测出来的基因序列进行聚类，构建非冗余基因集，得到非冗余基因集基因的碱基序列。

7. SOAPaligner 最大/最小插入片段长度：500/300 bp 基因丰度计算相似度（Identity）≥ 0.95。针对SOAPaligner比对后的信息，统计基因在各个样本中的丰度信息。

8. Diamond 比对类型： blastp E-value ≤ 1E-5，NR物种注释基于基因的物种分类学注释，比对NR数据库获得样本物种的分类学注释信息。

9. COG功能注释：Diamond 比对类型： blastp E-value ≤ 1E-5，比对EggNOG（evolutionary genealogy of genes: Non-supervised Orthologous Groups ）数据库获得基因对应的COG注释概况并进行统计

10. KEGG功能注释：Diamond 比对类型： blastp E-value ≤ 1E-5，比对KEGG（Kyoto Encyclopedia of Genes and Genomes）数据库获得基因对应的KEGG注释概况并进行统计。

11. CAZy碳水化合物活性酶注释，hmmscan 比对类型： hmmer E-value ≤ 1E-5，比对CAZy数据库（Carbohydrate-Active enZYmes Database）获得碳水化合物活性酶基因注释概况并进行统计。

12. ARDB抗性基因功能注释，Diamond 比对类型： blastp E-value ≤ 1E-5，比对ARDB（Antibiotic Resistance Genes Database）数据库获得抗性基因基因注释概况并进行统计。

13. CARD抗性基因功能注释，Diamond 比对类型： blastp E-value ≤ 1E-5，比对CARD（Comprehensive Antibiotic Resistance Database）数据库获得抗性基因注释概况并进行统计。

14. VFDB毒力因子注释，Diamond 比对类型： blastp E-value ≤ 1E-5，比对VFDB数据库获得毒力因子基因注释概况并进行统计。



http://picrust.github.io/picrust/tutorials/metagenome_prediction.html
