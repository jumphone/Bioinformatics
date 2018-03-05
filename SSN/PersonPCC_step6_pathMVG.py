import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.MVG','w')

ORIGIN=sys.argv[1].split('.')[-2].split('path2')[1] #Mtor
fa=open(sys.argv[1].split('.gene.path2Mtor.')[0]   ) #KO_sub_WT



INFO={}

for line in fa:
    seq=line.rstrip().split('\t')
    p1=seq[0]
    p2=seq[2]
   
    info='\t'.join(seq[4:8])
    if p1==ORIGIN:
        INFO[p2]=info
    elif p2==ORIGIN:
        INFO[p1]=info
    






GENE={}

for line in fi:
    seq=line.rstrip().split('\t')
    path=seq[2].split(',')
    if len(path)>2:
        #for gene in path[1:-1]:
            gene=path[1]
            if gene in GENE:
                GENE[gene]+=1
            else:
                GENE[gene]=1



output=[]

for gene in GENE:
    if gene in INFO:
        info=INFO[gene]
    else:
        info='NONE'
    output.append([GENE[gene],gene,info])

output.sort(reverse=True)

for one in output:
    fo.write(one[1]+'\t'+str(one[0])+'\t'+one[2]+'\n')
