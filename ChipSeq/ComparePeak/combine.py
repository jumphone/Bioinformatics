import sys
fa=open(sys.argv[1])
f1=open(sys.argv[2])
f2=open(sys.argv[3])
fo=open(sys.argv[4],'w')
old=[]
for line in fa:
    old.append(line.rstrip())


set1=set()
for line in f1:
    seq=line.rstrip().split('\t')
    if int(seq[-4])>0:
        set1.add(seq[0]+'\t'+seq[1]+'\t'+seq[2]+'\t'+seq[3])
        

set2=set()
for line in f2:
    seq=line.rstrip().split('\t')
    if int(seq[-4])>0:
        set2.add(seq[0]+'\t'+seq[1]+'\t'+seq[2]+'\t'+seq[3])

for one in old:
    f1t=0
    f2t=0
    if one in set1:
        f1t=1
    if one in set2:
        f2t=1


    fo.write(one+'\t'+str(f1t)+'\t'+str(f2t)+'\n')
