fi=open('10X_Tumor72017_matrix_CPTT.txt')
seq=fi.readline().rstrip().split('\t')
i=1
while i <len(seq):
    seq[i]='S10X.'+seq[i]
    i+=1
H1='\t'.join(seq[1:])
GENE=set()
G1={}
for line in fi:
    seq=line.replace('\n','').split('\t') 
    i=0
    while i<len(seq):
       if seq[i]=='':
            seq[i]='0.0'
       i+=1
    G1[seq[0]]='\t'.join(seq[1:])
    GENE.add(seq[0])
    


fi=open('GBM_normalized.txt')
seq=fi.readline().rstrip().split('\t')
i=1
while i <len(seq):
    seq[i]='SGBM.'+seq[i]
    i+=1
H2='\t'.join(seq[1:])
G2={}
for line in fi:
    seq=line.replace('\n','').split('\t')
    i=0
    while i<len(seq):
       if seq[i]=='':
            seq[i]='0.0'
       i+=1
    G2[seq[0]]='\t'.join(seq[1:])
    GENE.add(seq[0])







fo=open('COMBINED.txt','w')
H='GENE\t'+H1+'\t'+H2+'\n'
fo.write(H)
for G in GENE:
    fo.write(G+'\t')
    if G in G1 and G in G2:
        fo.write(G1[G])
        fo.write('\t')
        fo.write(G2[G])
        fo.write('\n')    














