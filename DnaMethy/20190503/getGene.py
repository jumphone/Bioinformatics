import  sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.txt','w')

GENE={}

for line in fi:
    seq=line.rstrip().split('\t')
    this_gene=seq[4]
    this_length=int(seq[2])-int(seq[1])
    if this_gene in GENE:
        GENE[this_gene] += float(seq[-1])/float(this_length)
    else:
        GENE[this_gene] = float(seq[-1])/float(this_length)


for gene in GENE:
    fo.write(gene+'\t'+str(GENE[gene])+'\n')

