import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
GL={}
GR={}
GP={}
TR=0
for line in fi:
    seq=line.rstrip().split('\t')
    this_gene=seq[4]
    this_l=float(seq[8])-float(seq[7])
    this_r=float(seq[-4])
    TR+=this_r

    if this_gene in GL and this_gene in GR:
        GP[this_gene]+=1
        GR[this_gene]+=this_r
        GL[this_gene]+=this_l
    else:
        GP[this_gene]=1
        GR[this_gene]=this_r
        GL[this_gene]=this_l

fo.write('gene\tpeaknumber\tlength\tread\tRPKM\n')
for gene in GL:
    rpkm = GR[gene] /( GL[gene] / 1000.0 * float(TR) /1000000.0  )
    fo.write(gene+'\t'+str(GP[gene])+'\t'+str(GL[gene])+'\t'+str(GR[gene])+'\t'+str(rpkm)+'\n')
