f1=open('CTRL.txt')
f2=open('KO.txt')
F2={}
for line in f2:
    seq=line.rstrip().split('\t')
    F2[seq[0]]='\t'.join(seq[1:])
fo=open('combined.txt','w')
old=set()
for line in f1:
    seq=line.rstrip().split('\t')
    if seq[0] in F2 and seq[0] not in old:
        old.add(seq[0])
        fo.write(line.rstrip()+'\t'+F2[seq[0]]+'\n')
