import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
i=1
for line in fi:
    seq=line.rstrip().split('\t')
    fo.write(seq[0]+'\t'+seq[1]+'\t'+seq[2]+'\tregion_'+str(i)+'\n')
    i=i+1
