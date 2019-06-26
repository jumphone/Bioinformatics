import sys
L=5000
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.5kb.bed','w')
for line in fi:
    seq=line.rstrip().split('\t')
    seq[1]=str(max(0,int(seq[1])-L))
    seq[2]=str(int(seq[2])+L)
    fo.write('\t'.join(seq)+'\n')
