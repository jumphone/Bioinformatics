'''
KRAS_1052=$wp\1052_S71_L007_R1_001.fastq.tophatdir/accepted_hits.bam
KRAS_1054=$wp\1054_S72_L007_R1_001.fastq.tophatdir/accepted_hits.bam
KRAS_1078=$wp\1078_S75_L007_R1_001.fastq.tophatdir/accepted_hits.bam
KRAS_1089=$wp\1089_S77_L007_R1_001.fastq.tophatdir/accepted_hits.bam
CONTROL_1055=$wp\1055_S73_L007_R1_001.fastq.tophatdir/accepted_hits.bam
CONTROL_1061=$wp\1061_S74_L007_R1_001.fastq.tophatdir/accepted_hits.bam
CONTROL_1081=$wp\1081_S76_L007_R1_001.fastq.tophatdir/accepted_hits.bam
CONTROL_1090=$wp\1090_S78_L007_R1_001.fastq.tophatdir/accepted_hits.bam
'''

KRAS_1052='1052_S71_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'
KRAS_1054='1054_S72_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'
KRAS_1078='1078_S75_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'
KRAS_1089='1089_S77_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'
CONTROL_1055='1055_S73_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'
CONTROL_1061='1061_S74_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'
CONTROL_1081='1081_S76_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'
CONTROL_1090='1090_S78_L007_R1_001.fastq.cufflinks/genes.fpkm_tracking'

FS=[KRAS_1052,KRAS_1054,KRAS_1078,KRAS_1089,CONTROL_1055,CONTROL_1061,CONTROL_1081,CONTROL_1090]
H=['GENE','KRAS_1052','KRAS_1054','KRAS_1078','KRAS_1089','CONTROL_1055','CONTROL_1061','CONTROL_1081','CONTROL_1090']
GS=[]
GENE=set()
for F in FS:
    G={}
    fi=open(F)
    fi.readline()
    for line in fi:
        seq=line.rstrip().split('\t')
        G[seq[4]]=seq[9]
        GENE.add(seq[4])
    GS.append(G)
fo=open('COMBINE.tsv','w')
header='\t'.join(H)+'\n'
fo.write(header)
for gene in GENE:
    tag=1
    for G in GS:
        if gene not in G:
            tag=0
    if tag==1:
        fo.write(gene)
        for G in GS:
            fo.write('\t'+G[gene])
        fo.write('\n')
