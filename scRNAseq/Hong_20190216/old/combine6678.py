fa=open('E-MTAB-6678.processed.1/E-MTAB-6687_arrayexpress_ss2_meta.txt')
header=fa.readline().rstrip().split('\t')
C2C={}
CC=set()
for line in fa:
    seq=line.rstrip().split('\t')
    C2C[seq[0]]=seq[2]
    CC.add(seq[2])

fi=open('E-MTAB-6678.processed.2/SS2_raw_counts.csv')
tmp=fi.readline()
G=set()
for line in fi:
    seq=line.rstrip().split(',')
    G.add(seq[0])


CLUSTER={}
for one in CC:
    CLUSTER[one]={}
    for g in G:
        CLUSTER[one][g]=0



fi=open('E-MTAB-6678.processed.2/SS2_raw_counts.csv')
header=fi.readline().rstrip().split(',')
for line in fi:
    seq=line.rstrip().split(',')
    this_gene=seq[0]
    i=1
    while i <len(seq):
        this_cell=header[i]
        if this_cell in C2C:
            this_cluster=C2C[this_cell]
            CLUSTER[this_cluster][this_gene]+=int(seq[i])
        i=i+1


H=[]
for one in CC:
    H.append(one)

R=[]
for one in G:
    R.append(one)

fo=open('ref_6678.txt','w')
fo.write('gene\t'+'\t'.join(H)+'\n')
for this_gene in R:
    fo.write(this_gene)
    for this_cluster in H:
        fo.write('\t'+str(CLUSTER[this_cluster][this_gene]))
    fo.write('\n')







