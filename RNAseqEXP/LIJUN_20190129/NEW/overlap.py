S2G={}
fi=open('whole_name.txt')
for line in fi:
    seq=line.rstrip().split('\t')
    S2G[seq[0]]=seq[1]

OPC=set()
NFOL=set()
MOL=set()
AC=set()
NEU=set()
MIC=set()

fi=open('ALL_TYPE.txt')
ref_header=fi.readline().rstrip()
for line in fi:
    seq=line.rstrip().split('\t')
    OPC.add(seq[0])
    NFOL.add(seq[1])
    MOL.add(seq[2])
    AC.add(seq[3])
    NEU.add(seq[4])
    MIC.add(seq[5])

fi=open('Brg1.txt')
header=fi.readline().rstrip()
fo=open('Brg1_combine.txt','w')

fo.write(header+'\t'+ref_header+'\twhole_name\n')

for line in fi:
    seq=line.rstrip().split('\t')
    sym=seq[1]
    wna='NA'
    if sym in S2G:
        wna=S2G[sym]

    opc='0'
    nfol='0'
    mol='0'
    ac='0'
    neu='0'
    mic='0'
    if sym in OPC:
        opc='1'
    if sym in NFOL:
        nfol='1'
    if sym in MOL:
        mol='1'
    if sym in AC:
    	ac='1'
    if sym in NEU:
    	neu='1'
    if sym in MIC:
    	mic='1'

    fo.write(line.rstrip()+'\t'+opc+'\t'+nfol+'\t'+mol+'\t'+ac+'\t'+neu+'\t'+mic+'\t'+wna+'\n')

