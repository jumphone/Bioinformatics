S2G={}
fi=open('while_name.txt')
for line in fi:
    seq=line.rstrip().split('\t')
    S2G[seq[0]]=seq[1]

OPC=set()
NFOL=set()
MOL=set()


fi=open('MARKER.txt')
fi.readline()
for line in fi:
    seq=line.rstrip().split('\t')
    OPC.add(seq[0])
    NFOL.add(seq[1])
    MOL.add(seq[2])


fi=open('Brg1 gene_exp.diff')
fo=open('Brg1 gene_exp.diff.combine.txt','w')

fo.write('Gene\tWhole\tSample1\tSample2\tExp1\tExp2\tPvalue\tQvalue\tOPC\tNFOL\tMOL\n')
fi.readline()
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[2] in S2G and seq[-1]=='yes':
        sym=seq[2]
        wna=S2G[seq[2]]
        sta1=seq[4]
        sta2=seq[5]
        v1=seq[7]
        v2=seq[8]
        pv=seq[11]
        qv=seq[12]
        opc='0'
        nfol='0'
        mol='0'
        if sym in OPC:
            opc='1'
        if sym in NFOL:
            nfol='1'
        if sym in MOL:
            mol='1'
        fo.write(sym+'\t'+wna+'\t'+sta1+'\t'+sta2+'\t'+v1+'\t'+v2+'\t'+pv+'\t'+qv+'\t'+opc+'\t'+nfol+'\t'+mol+'\n')

