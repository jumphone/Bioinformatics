import sys
anno_file='Homo_sapiens.GRCh37.75.chr.gtf'
input_file='REL_EXP.txt'
output_file='REL_EXP.sorted.txt'


fa=open(anno_file)
GENE={}
for line in fa:
    #if "protein_coding" in line and 'gene_id' in line and 'gene_name' in line and "chrN" not in line:
    if  'gene_id' in line and 'gene_name' in line and "chrN" not in line:
        seq=line.rstrip().split('\t')
        info=seq[8]
        gene_id=seq[8].split('gene_id "')[1].split('";')[0]
        gene_name=seq[8].split('gene_name "')[1].split('";')[0]

        if gene_id in GENE:
            GENE[gene_id][1].append(int(seq[3]))
            GENE[gene_id][1].append(int(seq[4]))
        else:
            GENE[gene_id]=[seq[0],[],gene_name]
            GENE[gene_id][1].append(int(seq[3]))
            GENE[gene_id][1].append(int(seq[4]))

for gene in GENE:
    GENE[gene][1].sort()

fi=open(input_file)
header=fi.readline().rstrip()
line=fi.readline().rstrip()
if line.count('\t') == header.count('\t')+1:
    header = 'gene\t'+header

output=[]
fi=open(input_file)
fi.readline()
for line in fi:
    seq=line.rstrip().split('\t')
    gene_id=seq[0]
    if gene_id in GENE:
        chrr=GENE[gene_id][0]
        start=GENE[gene_id][1][0]
        name=GENE[gene_id][2]
        output.append([chrr,start,name,gene_id+'\t'+'\t'.join(seq[1:])])

output.sort()
header='CHR\tSTART\tNAME\t'+header+'\n'
fo=open(output_file,'w')
fo.write(header)
for one in output:
    fo.write(one[0]+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\n')


   
