fa=open('E6701_2.csv')
header=fa.readline().rstrip().split('\t')
C2C={}
CC=set()
for line in fa:
    seq=line.rstrip().split('\t')
    C2C[seq[0]]=seq[1]+'\t'+seq[2]+'\t'+seq[3]
    CC.add(seq[2])


fo1=open('E6701_cluster.txt','w')



fi=open('E6701_1.txt')
header=fi.readline().rstrip().split('\t')
V=[]
i=1
while i<len(header):
    this_cell=header[i]
    if this_cell in C2C:
        fo1.write(this_cell+'\t'+C2C[this_cell]+'\n')
        V.append(i)
    i=i+1





fo2=open('E6701_exp.txt','w')
nh=['gene']
for i in V:
    nh.append(header[i])
fo2.write('\t'.join(nh)+'\n')


j=1
old=set()
for line in fi:
    if j %100==1:
        print(j)
    j=j+1
    seq=line.rstrip().split('\t')
    this_gene=seq[0].split('_')[0]
    if this_gene not in old:
        old.add(this_gene)
        out=[this_gene]
        for i in V:
            out.append(seq[i])
        fo2.write('\t'.join(out)+'\n')

