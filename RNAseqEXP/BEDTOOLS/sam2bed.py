import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[2]!='*':
        fo.write(seq[2]+'\t'+str(int(seq[3])-1)+'\t'+seq[3] +'\n')
fo.close()
