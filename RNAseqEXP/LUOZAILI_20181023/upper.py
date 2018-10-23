fi=open('NSC-KOCTD.txt')
fo=open('NSC-KOCTD.txt.tmp','w')
fo.write(fi.read().replace('\r','\n'))
fo.close()
fi=open('NSC-KOCTD.txt.tmp')
fo=open('NSC-KOCTD.txt.pure','w')
old=set()
for line in fi:
    seq=line.rstrip().split('\t')
    seq[0]=seq[0].upper()
    if seq[0] not in old:
        fo.write('\t'.join(seq)+'\n')
        old.add(seq[0])
