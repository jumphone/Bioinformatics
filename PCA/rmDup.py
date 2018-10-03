import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.rmdup','w')
old=set()
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[0] not in old:
        fo.write(line)
        old.add(seq[0])
