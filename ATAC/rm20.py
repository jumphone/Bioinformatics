fi=open('ALL_PEAK.sorted.merged.named.bed.p5pdRPKM.txt')
S=[]
for line in fi:
    seq=line.rstrip().split('\t')
    S.append(float(seq[6])+float(seq[7]))

S.sort()
C=S[int(len(S)*0.2)]
print(C)
fi=open('ALL_PEAK.sorted.merged.named.bed.kowtRPKM.txt')
fo=open('ANALYZED_RPKM.txt','w')
for line in fi:
    seq=line.rstrip().split('\t')
    if float(seq[6])+float(seq[7]) > C:
        fo.write(line)
