import sys
fi=open(sys.argv[1])
name=sys.argv[2]
fr=open(sys.argv[1]+'.barcode','w')
fo=open(sys.argv[1]+'.named.txt','w')

header=fi.readline().rstrip().split('\t')
i=1
while i<len(header):
    newname=name+'_'+str(i)
    fr.write(header[i]+'\t'+newname+'\n')
    header[i]=newname
    i+=1
fo.write('\t'.join(header)+'\n')
for line in fi:
    seq=line.rstrip().split('\t')
    i=1
    while i<len(seq):
        try:
            float(seq[i])
        except Exception as e:
            seq[i]=str(0.0)
        i+=1
    fo.write('\t'.join(seq)+'\n')



    





