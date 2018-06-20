import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.bedGraph','w')
fi.readline()
fo.write('track type=bedGraph\n')
for line in fi:
    seq=line.rstrip().split('\t')
    out=seq[1]+'\t'+str(int(seq[2])-1)+'\t'+seq[2]+'\t'+seq[5]+'\n'
    fo.write(out)
