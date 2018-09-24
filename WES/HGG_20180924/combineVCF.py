import sys

ALL=set()

SNV1={}
fi=open(sys.argv[1])
for line in fi:
    if line[0]!='#':
        seq=line.rstrip().split('\t')
        tag=seq[0]+':'+seq[1]+':'+seq[3]+':'+seq[4]
        info=line.rstrip()
        SNV1[tag]=info
        ALL.add(tag)


SNV2={}
fi=open(sys.argv[2])
for line in fi:
    if line[0]!='#':
        seq=line.rstrip().split('\t')
        tag=seq[0]+':'+seq[1]+':'+seq[3]+':'+seq[4]
        info=line.rstrip()
        SNV2[tag]=info
        ALL.add(tag)

SNV3={}
fi=open(sys.argv[3])
for line in fi:
    if line[0]!='#':
        seq=line.rstrip().split('\t')
        tag=seq[0]+':'+seq[1]+':'+seq[3]+':'+seq[4]
        info=line.rstrip()
        SNV3[tag]=info
        ALL.add(tag)

fo=open(sys.argv[4],'w')
OUT=[]
for tag in ALL:
    PASS=0
    if tag in SNV1:
        PASS+=1
        INFO=SNV1[tag]
    if tag in SNV2:
        PASS+=1
        INFO=SNV2[tag]
    if tag in SNV3:
        PASS+=1
        INFO=SNV3[tag]
    if PASS>=2:
        CHR=tag.split(':')[0]
        LOC=int(tag.split(':')[1])
        #OUT.append( [CHR,LOC,INFO+'\t'+str(PASS)])
        OUT.append( [CHR,LOC,INFO])
OUT.sort()
for one in OUT:
    fo.write(one[2]+'\n')
fo.close()
