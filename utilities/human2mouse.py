import sys
human2HGNC={}
HGNC2mouse={}



fa=open('./hgnc_complete_set.txt')
for line in fa:
    if "protein-coding" in line:

        seq=line.rstrip().split('\t')
        hgnc=seq[0]
        symbol=seq[1]
        human2HGNC[symbol]=hgnc



fa=open('./MGI_HGNC_homologene.rpt')
for line in fa:
    seq=line.rstrip().split('\t')
    mouse=seq[1]
    ens=seq[9]
    hgnc=seq[20]
    HGNC2mouse[hgnc]=mouse+'\t'+ens



#print human2HGNC


fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.mouse','w')
for line in fi:
    seq=line.rstrip().split('\t')
    human=seq[0].strip()
    #print human
    #if human in human2HGNC:
    #    print human
    if human in human2HGNC and human2HGNC[human] in HGNC2mouse:
        hgnc=human2HGNC[human]
        mouse=HGNC2mouse[hgnc]
        fo.write(human+'\t'+hgnc+'\t'+mouse+'\n')
