
python 1changeName.py ../data_new/1.5M_D_out_gene_exon_tagged.dge.txt S1.5MD 
python 1changeName.py ../data_new/1.5M_P_out_gene_exon_tagged.dge.txt S1.5MP
python 1changeName.py ../data_new/1.5M_SN_out_gene_exon_tagged.dge.txt S1.5MSN
python 1changeName.py ../data_new/4M_D_out_gene_exon_tagged.dge.txt S4MD
python 1changeName.py ../data_new/4M_P_out_gene_exon_tagged.dge.txt S4MP
python 1changeName.py ../data_new/4M_SN_out_gene_exon_tagged.dge.txt S4MSN


DATA='../data_new/'
python 2combineMatrix.py $DATA\1.5M_D_out_gene_exon_tagged.dge.txt.named.txt $DATA\1.5M_P_out_gene_exon_tagged.dge.txt.named.txt $DATA\EXP12.combined.txt
python 2combineMatrix.py $DATA\EXP12.combined.txt $DATA\1.5M_SN_out_gene_exon_tagged.dge.txt.named.txt $DATA\EXP123.combined.txt
python 2combineMatrix.py $DATA\EXP123.combined.txt $DATA\4M_D_out_gene_exon_tagged.dge.txt.named.txt $DATA\EXP1234.combined.txt
python 2combineMatrix.py $DATA\EXP1234.combined.txt $DATA\4M_P_out_gene_exon_tagged.dge.txt.named.txt $DATA\EXP12345.combined.txt
python 2combineMatrix.py $DATA\EXP12345.combined.txt $DATA\4M_SN_out_gene_exon_tagged.dge.txt.named.txt $DATA\EXP123456.combined.txt


#python 3cellExpGene.py ../data/EXP123456.combined.txt
#python 4varEXP.py ../data/EXP123456.combined.txt 

