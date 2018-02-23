import sys
f1=open(sys.argv[1])
f2=open(sys.argv[2])
fo=open(sys.argv[3],'w')
GENE={}

header1=f1.readline().rstrip().split('\t')[1:]
header2=f2.readline().rstrip().split('\t')[1:]
newheader='GENE\t'+'\t'.join(header1)+'\t'+'\t'.join(header2)+'\n'
fo.write(newheader)

GENE1={}
for line in f1:
    seq=line.rstrip().split('\t')
    GENE1[seq[0]]=seq[1:]
    GENE[seq[0]]=['0']*(len(header1)+len(header2))
GENE2={}
for line in f2:
    seq=line.rstrip().split('\t')
    GENE2[seq[0]]=seq[1:]
    GENE[seq[0]]=['0']*(len(header1)+len(header2))

for gene in GENE1:
    GENE[gene][0:len(header1)]= GENE1[gene]

for gene in GENE2:
    GENE[gene][len(header1):]= GENE2[gene]
