
i2g={}
fi=open('INTER5kb.txt')
for line in fi:
    seq=line.rstrip().split('\t')
    this_i='.'.join(seq[:3])
    #print(this_i)
    if this_i in i2g:
        if seq[6] not in i2g[this_i]:
            i2g[this_i].append(seq[6])
    else:
         i2g[this_i]=[seq[6]]



INDEX=[]
RPKM={}
fi=open('Lijun_ATAC_1.rmdup.bam.bed.cov.rpkm')
for line in fi:
    seq=line.rstrip().split('\t')
    this_i='.'.join(seq[:3])
    INDEX.append(this_i)
    RPKM[this_i]=[seq[7]]
    


fi=open('Lijun_ATAC_2.rmdup.bam.bed.cov.rpkm')
for line in fi:
    seq=line.rstrip().split('\t')
    this_i='.'.join(seq[:3])
    RPKM[this_i].append(seq[7])


fi=open('Lijun_ATAC_3.rmdup.bam.bed.cov.rpkm')
for line in fi:
    seq=line.rstrip().split('\t')
    this_i='.'.join(seq[:3])
    RPKM[this_i].append(seq[7])

fi=open('Lijun_ATAC_4.rmdup.bam.bed.cov.rpkm')
for line in fi:
    seq=line.rstrip().split('\t')
    this_i='.'.join(seq[:3])
    RPKM[this_i].append(seq[7])


header=['CHR','START','END','LIJUN1','LIJUN2','LIJUN3','LIJUN4','GENE']
fo=open('COMBINE_5kb.txt','w')
fo.write('\t'.join(header)+'\n')
for one in INDEX:
    gene=''
    if one in i2g:
        gene=','.join(i2g[one])
    fo.write(one.replace('.','\t')+'\t'+'\t'.join(RPKM[one])+'\t'+  gene+'\n' )








