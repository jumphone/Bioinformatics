import sys
HGNC2human={}
mouse2HGNC={}

fa=open('./hgnc_complete_set.txt')
for line in fa:
    if "protein-coding" in line:

        seq=line.split('\t')
        hgnc=seq[0]
        symbol=seq[1]
        HGNC2human[hgnc]=symbol



fa=open('./MGI_HGNC_homologene.rpt')
for line in fa:
    seq=line.rstrip().split('\t')
    sym=seq[1]
    ens=seq[9]
    hgnc=seq[20]
    mouse2HGNC[sym]=hgnc
    

fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.human','w')
for line in fi:
    seq=line.rstrip().split('\t')
    mouse=seq[0]
    if mouse in mouse2HGNC and mouse2HGNC[mouse] in HGNC2human:
        hgnc=mouse2HGNC[mouse]
        human=HGNC2human[hgnc]
        fo.write(mouse+'\t'+hgnc+'\t'+human+'\t'+line.rstrip()+'\n')

    
