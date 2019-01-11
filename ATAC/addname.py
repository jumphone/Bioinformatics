fi=open('ALL_PEAK.sorted.merged.bed')
fo=open('ALL_PEAK.sorted.merged.named.bed','w')

i=1
for line in fi:
    NAME='ALLPEAK_'+str(i)
    fo.write(line.rstrip()+'\t'+NAME+'\n')
    i+=1
