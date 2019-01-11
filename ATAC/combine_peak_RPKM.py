import sys
P5=sys.argv[1]#'P5-50K-jincheng.rmdup.bam.ALLPEAK'
PD=sys.argv[2]#'PDGFRa-GFP-50K-E14.5-jincheng.rmdup.bam.ALLPEAK'

P5_reads_in_peaks=0
PD_reads_in_peaks=0


fi=open(P5)
P5_reads={}
for line in fi:
    seq=line.rstrip().split('\t')
    P5_reads_in_peaks+=int(seq[-4])
    P5_reads[seq[-5]]=int(seq[-4])


fi=open(PD)
PD_reads={}
for line in fi:
    seq=line.rstrip().split('\t')
    PD_reads_in_peaks+=int(seq[-4])
    PD_reads[seq[-5]]=int(seq[-4])



fa=open('Mus_musculus.GRCm38.87.chr.gtf.combined.pc.bed.sort.tss.bed')
TSS={}
for line in fa:
    seq=line.rstrip().split('\t')
    tag=seq[3]+'\t'+seq[4]
    chrr=seq[0]
    loc=int(seq[2])
    try:
        TSS[chrr].append([loc,tag])
    except Exception as e:
        TSS[chrr]=[[loc,tag]]




print P5_reads_in_peaks
print PD_reads_in_peaks

#P5_N=206579
#PD_N=1115469
#print len(P5_reads)
#print len(PD_reads)


fi=open('ALL_PEAK.sorted.merged.named.bed')
fo=open('ALL_PEAK.sorted.merged.named.bed.kowtRPKM.txt','w')
i=1
for line in fi:
    print i
    i+=1
    seq=line.rstrip().split('\t')
    #print seq
    p5=P5_reads[seq[3]]
    pd=PD_reads[seq[3]]
    CHRR=seq[0]
    CENTER=int((int(seq[2])+int(seq[1])+1.0)/2.0)
    lll=float(int(seq[2])-int(seq[1]))
    DIS=[]
    for tss in TSS[CHRR]:
        DIS.append([abs(tss[0]-CENTER),tss[0],tss[1]])
    DIS.sort()


    P5_rpkm=(float(p5))/(float(P5_reads_in_peaks)/1000000.0*lll/1000.0)+1
    PD_rpkm=(float(pd))/(float(PD_reads_in_peaks)/1000000.0*lll/1000.0)+1
    fo.write(line.rstrip()+'\t'+str(p5)+'\t'+str(pd)+'\t'+str(P5_rpkm)+'\t'+str(PD_rpkm)+'\t'+str(P5_rpkm/PD_rpkm)+'\t'+str(DIS[0][0])+'\t'+DIS[0][2]+'\t'+str(DIS[0][1])+'\n')


