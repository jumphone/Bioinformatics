import cPickle as pickle
import sys
from scipy import stats

import multiprocessing

print '''
$1: Index Data
$2: EXP Matrix
$3: Output file
$4: CPU
'''

INDEX_DATA_FILE=sys.argv[1]
EXP_MATRIX=sys.argv[2]
OUT_FILE=sys.argv[3]
PROC_LIMIT=int(sys.argv[4])



class Person_PCC_Data:
    def __init__(self,  EXP , POOL_LENGTH , PCC_POOL ):
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL

print 'loading...'
fdata = open(INDEX_DATA_FILE)
data = pickle.load(fdata)
fdata.close()
print "loading done !"


THIS_EXP={}
fi=open(EXP_MATRIX)
line=fi.readline()
header=line.rstrip().split('\t')

for one in header:
    THIS_EXP[one]={}

for line in fi:
    try:
        seq=line.rstrip().split('\t')
        i=1
        while i<len(seq):
            H=header[i-1]
            THIS_EXP[H][seq[0]]=float(seq[i])
            i+=1
    except Exception as e:
        pass
fi.close()
import random,subprocess
OUTDIR=OUT_FILE+'_TMP_'+str(random.random())
subprocess.Popen('mkdir '+OUTDIR,shell=True).wait()


def SINGLE(p1, p2,p1_old_exp,p2_old_exp):
    outfile=OUTDIR+'/'+p1+'_'+p2+'.ssn'
    if os.path.exists(outfile)!=True:
        Z=[]
        j=0
        while j<len(header):
            H=header[j]
            p1_new_exp=p1_old_exp+[Person_EXP[H][p1]]
            p2_new_exp=p2_old_exp+[Person_EXP[H][p2]]
            pcc_new = stats.pearsonr(p1_new_exp,p2_new_exp)[0]
            pcc = data.PCC_POOL[edge]
            delta_pcc = pcc_new - pcc
            z=delta_pcc /( (1-pcc**2)/(data.POOL_LENGTH-1) )
            Z.append(str(z))
            j+=1
        fo=open(outfile,'w')
        fo.write(p1+'.And.'+p2+'\t'+'\t'.join(Z)+'\n')



open(OUTDIR+'header.txt','w').write('\t'.join(header)+'\n')

jobs=[]
for edge in data.PCC_POOL:
    ps=edge.split(':')
    p1=ps[0]
    p2=ps[1]
    try:
        p1_old_exp=data.EXP[p1]
        p2_old_exp=data.EXP[p2]
        Person_EXP[header[0]][p1]
        Person_EXP[header[0]][p2]
        p= multiprocessing.Process(target=SINGLE, args=(p1,p2,p1_old_exp,p2_old_exp))
        p.start()
        jobs.append(p)
    except Exception as e:
        pass
    if len(jobs)>=PROC_LIMIT:
        for p in jobs:
            p.join()
        jobs=[]
for p in jobs:
    p.join()

subprocess.Popen('cat '+OUTDIR+'/*.ssn > '+OUTDIR+'/DATA.SSN',shell=True).wait()
subprocess.Popen('cat '+OUTDIR+'/header.txt '+OUTDIR+'/DATA.SSN > '+OUT_FILE,shell=True).wait()
























