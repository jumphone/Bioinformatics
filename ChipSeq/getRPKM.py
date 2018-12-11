fwt=open('ALL.peak.sort.merged.WT.cov') #WT
fko=open('ALL.peak.sort.merged.KO.cov') #KO
fo=open('WT_KO_RPKM.txt','w')

fo.write('#CHR\tSTART\tEND\tWT\tKO\n')

WT=fwt.read().split('\n')
KO=fko.read().split('\n')
all_wt = 2553050 / 1000000.0
all_ko = 3639947 / 1000000.0

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


        output.append([wt[0],int(wt[1]),sss ,wt[0]+'\t'+wt[1]+'\t'+wt[2]+'\t'+str(rwt)+'\t'+str(rko)+'\n'])
    i+=1

output.sort()
for one in output:
    fo.write(one[3].rstrip()+'\n')
