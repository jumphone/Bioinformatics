import sys

fa=open('Homo_sapiens.GRCh37.75.chr.gtf')
name2id={}
for line in fa:
    if line[0]!='#' and 'gene_name' in line and 'gene_id' in line:
        name= line.split('gene_name "')[1].split('";')[0]
        gid= line.split('gene_id "')[1].split('";')[0]

        name2id[name]=gid

fi=open('GSE70630_OG_processed_data_v2.txt')
fo=open('GSE70630_OG_processed_data_v2.txt.pure','w')
fo.write(fi.readline())
for line in fi:
    seq=line.rstrip().split('\t')
    name=seq[0].replace("'",'')
    if name in name2id:
        gid= name2id[name]
        fo.write(gid+'\t'+'\t'.join(seq[1:])+'\n')
