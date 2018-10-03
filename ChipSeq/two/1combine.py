fwt=open('ALL.peak.sort.merged.WT.cov') #WT
fko=open('ALL.peak.sort.merged.KO.cov') #KO
fo=open('SEQ.WT_KO','w')

fo.write('#CHR\tSTART\tEND\tWT\tKO\tKO_div_WT\tSTART_NEW\tEND_NEW\n')

KO=fwt.read().split('\n')
WT=fko.read().split('\n')


all_wt = 3141421 / 1000000.0
all_ko = 2901128 / 1000000.0





i=0
output=[]
SSS=[]

def change(ratio):
    if ratio ==1:
        ratio = 0
    elif ratio <1:
        ratio = -1/ratio
    return ratio

while i<len(KO):
    ko=KO[i].split('\t')

    wt=WT[i].split('\t')
    if len(ko)>3:
        l= float(int(wt[2])-int(wt[1]))/1000.0

        rwt=float(wt[3])/l / all_wt
        rko=float(ko[3])/l / all_ko

        sss=sum([rwt,rko])
        SSS.append(sss)

        ko_div_wt = (rko+1) / (rwt+1)
        ko_div_wt = change(ko_div_wt)


        output.append([wt[0],int(wt[1]),sss ,wt[0]+'\t'+wt[1]+'\t'+wt[2]+'\t'+str(rwt)+'\t'+str(rko)+'\t'+str(ko_div_wt)+'\n'])
    i+=1

SSS.sort()
LIMIT=0.1
CUTOFF=SSS[int(LIMIT*len(SSS))]
output.sort()
print CUTOFF

tmp_chrr=""
tmp_end=0
#fo_wt=open('ALL.ATAC.narrowPeak.merged.WT_KO.WT','w')
#fo_ko=open('ALL.ATAC.narrowPeak.merged.WT_KO.KO','w')
KO_WT_u_KO={}
KO_WT_u_WT={}




CHRR_LIST=[]
CHRR_SUM={}
for one in output:
    if one[2] > CUTOFF:
        seq=one[3].rstrip().split('\t')
        chrr=seq[0]
        if chrr!='chrM' and chrr!='chrY':
            start = int(seq[1])
            end = int(seq[2])
            L=end-start
            if chrr !=tmp_chrr:
               CHRR_LIST.append(chrr)
               start_new=0
               end_new = start_new +L
            else:
               start_new = tmp_end
               end_new = start_new+L
            if chrr in CHRR_SUM:
                CHRR_SUM[chrr]=max([CHRR_SUM[chrr],end_new])
            else:
                CHRR_SUM[chrr]=end_new

            tmp_chrr = chrr
            tmp_end = end_new
            ko_div_wt=float(seq[5])



            if chrr in KO_WT_u_KO:
                if ko_div_wt > 1:
                    KO_WT_u_KO[chrr]+=L
            else:
                if ko_div_wt > 1:
                    KO_WT_u_KO[chrr]=L
            if chrr in KO_WT_u_WT:
                if ko_div_wt < -1:
                    KO_WT_u_WT[chrr]+=L
            else:
                if ko_div_wt < -1:
                    KO_WT_u_WT[chrr]=L

            fo.write(one[3].rstrip()+'\t'+str(start_new)+'\t'+str(end_new)+'\n')

fo.close()

fo1=open('SEQ.KO_WT_u_KO','w')
fo4=open('SEQ.KO_WT_u_WT','w')



for chrr in CHRR_LIST:
    fo1.write(chrr+'\t0\t'+str(KO_WT_u_KO[chrr])+'\t1\t'+str(int(KO_WT_u_KO[chrr]/float(CHRR_SUM[chrr])*100))+'\n')
    #print KO_WT_u_KO[chrr]
    #print KO_WT_u_WT[chrr]

    fo4.write(chrr+'\t'+str(KO_WT_u_KO[chrr])+'\t'+str(KO_WT_u_KO[chrr]+KO_WT_u_WT[chrr])+'\t0\t'+str(int(KO_WT_u_WT[chrr]/float(CHRR_SUM[chrr])*100))+'\n')






