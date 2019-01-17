import sys
P5=sys.argv[1]

fi=open(P5)
fo=open(sys.argv[2],'w')


P5_reads_in_peaks=0

P5_reads={}
for line in fi:
    seq=line.rstrip().split('\t')
    P5_reads_in_peaks+=int(seq[-4])




fi=open(P5)
for line in fi:
    seq=line.rstrip().split('\t')
    lll=float(seq[2])-float(seq[1])
    rpkm=(float(seq[-4]))/(float(P5_reads_in_peaks)/1000000.0*lll/1000.0)
    fo.write(seq[3]+'\t'+str(rpkm)+'\n')
