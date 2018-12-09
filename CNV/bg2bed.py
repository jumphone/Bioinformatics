import sys

fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.bed','w')

BED={}
CHR=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
        'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18',
        'chr19','chr20','chr21','chrX','chrY']
for chrr in CHR:
    BED[chrr]=[]


FLAG='NA'
for line in fi:
    if "track" in line:
        if "Gain" in line:
            FLAG='1'
        elif 'Loss' in line:
            FLAG='-1'
        else:
            FLAG='0'
    else:
        seq=line.rstrip().split('\t')
        if seq[0] in BED:

            BED[seq[0]].append([int(seq[1]),int(seq[2]),seq[3],FLAG])

for chrr in BED:
    BED[chrr].sort()
    for one in BED[chrr]:
        fo.write(chrr+'\t'+str(one[0])+'\t'+str(one[1])+'\t'+one[2]+'\t'+one[3]+'\n')


