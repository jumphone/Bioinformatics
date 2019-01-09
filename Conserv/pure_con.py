import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.pure','w')
for line in fi:
        seq=line.replace('\n','').split('\t')
        for one in seq[2:]:
                try:
                        float(one)
                        fo.write(one+'\n')
                except Exception:
                        next
