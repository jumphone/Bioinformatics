import gzip
fi=open('sgRNA.txt')

def antisense_reverse(read):
        read=read.upper()
        read_change_base=""
        for one in read:
                if one == 'A':
                        read_change_base += 'T'
                elif one == 'C':
                        read_change_base += 'G'
                elif one == 'G':
                        read_change_base += 'C'
                elif one == 'T':
                        read_change_base += 'A'
                else:
                        read_change_base += 'N'
        read_reverse=read_change_base[::-1]
        return read_reverse

SEQ={}
for line in fi:
    seq=line.rstrip().split('\t')
    this_seq=seq[3][1:-1]
    this_seq_rev=antisense_reverse(this_seq)
    SEQ[seq[0]]=[this_seq,this_seq_rev]


STAT={}
for sg in SEQ:
   STAT[sg]=0



import sys
fi=gzip.open(sys.argv[1],'rb')
fo=open(sys.argv[1]+'.STAT.txt','w')
l1=fi.readline()
l2=fi.readline()
l3=fi.readline()
l4=fi.readline()

i=1
while l1 !='':

    for sg in SEQ:
        if SEQ[sg][0] in l2 or SEQ[sg][1] in l2:
            STAT[sg] +=1
            #break
    if i%100==1:
        print(i)
    #if i >3000: 
    #    break
    i=i+1
    l1=fi.readline()
    l2=fi.readline()
    l3=fi.readline()
    l4=fi.readline()

for sg in STAT:
    fo.write(sg+'\t'+str(STAT[sg])+'\n')
