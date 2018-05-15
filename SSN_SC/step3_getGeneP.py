import sys
from scipy import stats
import multiprocessing
print '''
ZMAT=sys.argv[1]
OUT=sys.argv[2]
CPU=int(sys.argv[3])
'''

ZMAT=sys.argv[1]
OUT=sys.argv[2]
CPU=int(sys.argv[3])
#REFG=sys.argv[2]

fi=open(ZMAT)
header=fi.readline().rstrip().split('\t')
CELL={}
for cell in header:
    CELL[cell]={}
j=1
for line in fi:
    #print j;j+=1
    seq=line.rstrip().split('\t')
    g1=seq[0].split('.And.')[0]
    g2=seq[0].split('.And.')[1]
    i=1
    while i<len(seq):
        cell=header[i-1]
        try:
            zvalue=float(seq[i])
        except Exception as e:
            zvalue=0.0

        score= abs(zvalue)#stats.norm.sf(abs(zvalue))*2 #abs(zvalue)

        if 'ALL' in CELL[cell]:
            CELL[cell]['ALL'].append(score)
        else:
            CELL[cell]['ALL']=[score]
        if g1 in CELL[cell]:
            CELL[cell][g1].append(score)
        else:
            CELL[cell][g1]=[score]
        if g2 in CELL[cell]:
            CELL[cell][g2].append(score)
        else:
            CELL[cell][g2]=[score]
        i+=1


#for cell in CELL:

GENE_INDEX=[]
for gene in CELL[header[0]]:
    if gene !='ALL':
        GENE_INDEX.append(gene)
#print gene

import random,subprocess
OUTDIR=OUT+'_TMP_'+str(random.random())
subprocess.Popen('mkdir '+OUTDIR,shell=True).wait()



def SINGLE(tmp,gene):
    outfile=OUTDIR+'/'+gene+'.pvalue'
    fo=open(outfile,'w')
    fo.write(gene)
    for cell in header:
        try:
            pvalue= stats.mannwhitneyu( CELL[cell][gene], CELL[cell]['ALL'],use_continuity=True,alternative='greater')[1]
            #stats.kstest(CELL[cell][gene], 'norm')[1]
        except Exception as e:
            pvalue=1.0
        fo.write('\t'+str(pvalue))
    fo.write('\n')


fo=open(OUTDIR+'/header.txt','w')
fo.write('\t'.join(header)+'\n')
PROC_LIMIT=CPU
jobs=[]
for gene in GENE_INDEX:
    try:
        p= multiprocessing.Process(target=SINGLE, args=(1,gene))
        p.start()
        jobs.append(p)
    except Exception as e:
        print e
    if len(jobs)>=PROC_LIMIT:
        for p in jobs:
            p.join()
        jobs=[]
for p in jobs:
    p.join()

subprocess.Popen('cat '+OUTDIR+'/*.pvalue > '+OUTDIR+'/DATA.pvalue',shell=True).wait()
subprocess.Popen('cat '+OUTDIR+'/header.txt '+OUTDIR+'/DATA.pvalue > '+OUT,shell=True).wait()


#stats.mannwhitneyu([1,2,3],[3,4,5,7],use_continuity=True,alternative='less') #greater








