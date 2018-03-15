f1=open('ALL.narrowPeak.sort.merged.KO.cov')
f2=open('ALL.narrowPeak.sort.merged.KI.cov')
f3=open('ALL.narrowPeak.sort.merged.WT.cov')
fo=open('CHIP.WT_KO_KI','w')

fo.write('#CHR\tSTART\tEND\tWT\tKO\tKI\tKO_div_WT\tKI_div_WT\tKO_div_KI\tSTART_NEW\tEND_NEW\n')

KO=f1.read().split('\n')
KI=f2.read().split('\n')
WT=f3.read().split('\n')


all_wt = 26540355 / 1000000.0
all_ko = 25564230 / 1000000.0
all_ki = 23913357  / 1000000.0




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
    ki=KI[i].split('\t')
    wt=WT[i].split('\t')
    if len(ko)>3:
        l= float(int(wt[2])-int(wt[1]))/1000.0

        rwt=float(wt[3])/l / all_wt
        rko=float(ko[3])/l / all_ko
        rki=float(ki[3])/l / all_ki
        sss=sum([rwt,rko,rki])
        SSS.append(sss)

        ko_div_wt = (rko+1) / (rwt+1)
        ki_div_wt = (rki+1) / (rwt+1)
        ko_div_ki = (rko+1) / (rki+1)

        ko_div_wt = change(ko_div_wt)
        ki_div_wt = change(ki_div_wt)
        ko_div_ki = change(ko_div_ki)
        
        output.append([wt[0],int(wt[1]),sss ,wt[0]+'\t'+wt[1]+'\t'+wt[2]+'\t'+str(rwt)+'\t'+str(rko)+'\t'+str(rki)+'\t'+str(ko_div_wt)+'\t'+str(ki_div_wt)+'\t'+str(ko_div_ki)+'\n'])
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
KI_WT_u_KI={}
KI_WT_u_WT={}
KO_KI_u_KO={}
KO_KI_u_KI={}

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
            ko_div_wt=float(seq[6])
            ki_div_wt=float(seq[7])
            ko_div_ki=float(seq[8])


            if chrr in KO_WT_u_KO:
                if ko_div_wt > 1:
                    KO_WT_u_KO[chrr]+=L
            else:
                if ko_div_wt > 1:
                    KO_WT_u_KO[chrr]=L
            if chrr in KI_WT_u_KI:
                if ki_div_wt  > 1:
                    KI_WT_u_KI[chrr]+=L
            else:
                if ki_div_wt  > 1:
                    KI_WT_u_KI[chrr]=L
            if chrr in KO_KI_u_KO:
                if ko_div_ki > 1:
                    KO_KI_u_KO[chrr]+=L
            else:
                if ko_div_ki > 1:
                    KO_KI_u_KO[chrr]=L

            if chrr in KO_WT_u_WT:
                if ko_div_wt < -1:
                    KO_WT_u_WT[chrr]+=L
            else:
                if ko_div_wt < -1:
                    KO_WT_u_WT[chrr]=L
            if chrr in KI_WT_u_WT:
                if ki_div_wt  < -1:
                    KI_WT_u_WT[chrr]+=L
            else:
                if ki_div_wt  < -1:
                    KI_WT_u_WT[chrr]=L
            if chrr in KO_KI_u_KI:
                if ko_div_ki < -1:
                    KO_KI_u_KI[chrr]+=L
            else:
                if ko_div_ki < -1:
                    KO_KI_u_KI[chrr]=L

            fo.write(one[3].rstrip()+'\t'+str(start_new)+'\t'+str(end_new)+'\n')

fo.close()

fo1=open('CHIP.KO_WT_u_KO','w')
fo2=open('CHIP.KI_WT_u_KI','w')
fo3=open('CHIP.KO_KI_u_KO','w')
fo4=open('CHIP.KO_WT_u_WT','w')
fo5=open('CHIP.KI_WT_u_WT','w')
fo6=open('CHIP.KO_KI_u_KI','w')

for chrr in CHRR_LIST:
    fo1.write(chrr+'\t0\t'+str(KO_WT_u_KO[chrr])+'\t1\t'+str(int(KO_WT_u_KO[chrr]/float(CHRR_SUM[chrr])*100))+'\n')
    fo2.write(chrr+'\t0\t'+str(KI_WT_u_KI[chrr])+'\t1\t'+str(int(KI_WT_u_KI[chrr]/float(CHRR_SUM[chrr])*100))+'\n')
    fo3.write(chrr+'\t0\t'+str(KO_KI_u_KO[chrr])+'\t1\t'+str(int(KO_KI_u_KO[chrr]/float(CHRR_SUM[chrr])*100))+'\n')
    fo4.write(chrr+'\t'+str(KO_WT_u_KO[chrr])+'\t'+str(KO_WT_u_KO[chrr]+KO_WT_u_WT[chrr])+'\t0\t'+str(int(KO_WT_u_WT[chrr]/float(CHRR_SUM[chrr])*100))+'\n')
    fo5.write(chrr+'\t'+str(KI_WT_u_KI[chrr])+'\t'+str(KI_WT_u_KI[chrr]+KI_WT_u_WT[chrr])+'\t0\t'+str(int(KI_WT_u_WT[chrr]/float(CHRR_SUM[chrr])*100))+'\n')
    fo6.write(chrr+'\t'+str(KO_KI_u_KO[chrr])+'\t'+str(KO_KI_u_KO[chrr]+KO_KI_u_KI[chrr])+'\t0\t'+str(int(KO_KI_u_KI[chrr]/float(CHRR_SUM[chrr])*100))+'\n')











        
