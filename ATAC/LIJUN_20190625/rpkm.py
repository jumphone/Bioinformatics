import sys
fi=open(sys.argv[1])

fo=open(sys.argv[1]+'.rpkm','w')

NUM=0
for line in fi:
    seq=line.rstrip().split('\t')
    this_num=float(seq[-4])
    #print(this_num)
    NUM=NUM+this_num

fi.close()
fi=open(sys.argv[1])

for line in fi:
    seq=line.rstrip().split('\t')
    this_num=float(seq[-4])
    numReads=this_num
    geneLength=float(seq[2])-float(seq[1])
    totalNumReads=NUM
    rpkm=numReads / ( geneLength/1000 * totalNumReads/1000000 )
    rpkm=round(rpkm,2)
    fo.write('\t'.join(seq+[str(rpkm)])+'\n')



