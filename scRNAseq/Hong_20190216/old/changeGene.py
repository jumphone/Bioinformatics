G2G={}

fi=open('ref_6701.txt')
fo=open('ref_6701_pure.txt','w')
fo.write(fi.readline())
old=set()
for line in fi:
    seq=line.rstrip().split('\t')
    info=seq[0].split('_')
    this_gene=info[0]
    this_ens=info[1]
    G2G[this_ens]=this_gene
    if this_gene not in old:
        old.add(this_gene)
        seq[0]=this_gene
        fo.write('\t'.join(seq)+'\n')
fo.close()

fi=open('ref_6678.txt')
fo=open('ref_6678_pure.txt','w')
fo.write(fi.readline())
old=set()
for line in fi:
    seq=line.rstrip().split('\t')
    this_ens=seq[0]
    if this_ens in G2G:
        this_gene=G2G[this_ens]
        if this_gene not in old:
            old.add(this_gene)
            seq[0]=this_gene
            fo.write('\t'.join(seq)+'\n')
fo.close()
