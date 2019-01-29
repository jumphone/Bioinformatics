import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.rmNa','w')
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[0] !='.':
        fo.write(line)
