import sys
print '''
file1: step3 output

'''

fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.gene','w')

PAIR={}

for line in fi:
    seq=line.split('\t')
    p1=seq[0]
    p2=seq[2]
    tag=[p1,p2]
    tag.sort()
    tag=':'.join(tag)
    PAIR[tag]=float(seq[4])

GENE={}

for pair in PAIR:
    p1=pair.split(':')[0]
    p2=pair.split(':')[1]
    z=abs(PAIR[pair])
    if p1 in GENE:
        GENE[p1] += z
    else:
        GENE[p1]=z
    if p2 in GENE:
        GENE[p2] += z
    else:
        GENE[p2]=z

output=[]
for one in GENE:
    output.append([GENE[one],one])
i=1

output.sort(reverse=True)
for one in output:
    fo.write(str(i)+'\t'+one[1]+'\t'+str(one[0])+'\n')
    i+=1
