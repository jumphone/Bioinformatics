import sys
P5={}
E14={}
DISTANCE=5000
fi=open('ANALYZED_RPKM.txt')
for line in fi:
    seq=line.rstrip().split('\t')
    if int(seq[-4])<=DISTANCE:
        gene=seq[-2]
        p5=float(seq[6])
        e14=float(seq[7])
        if gene in P5:
            P5[gene]=P5[gene]+p5
            E14[gene]=E14[gene]+e14
        else:
            P5[gene]=p5
            E14[gene]=e14
fo=open('ANALYZED_RPKM.txt.ok.txt','w')
#OLD=set()
#faa=open('MGI.txt')
#for line in faa:
#    OLD.add(line.rstrip())

fo.write('NAME\tDESCRIPTION\tKO\tWT\n')
for gene in P5:
    #if gene in OLD:
        fo.write(gene+'\tna\t'+str(P5[gene])+'\t'+str(E14[gene])+'\n')
