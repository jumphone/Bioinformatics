import sys
import cPickle as pickle

print '''

$1 NETWORK_PATH (col1: gene; col2: gene)
$2 EXP_POOL_PATH (col 1: gene; other col: exp value)
$3 OUT_DATA_PATH

'''



NETWORK_PATH=sys.argv[1]
EXP_POOL_PATH=sys.argv[2]
OUT_DATA_PATH=sys.argv[3]

class Person_PCC_Data:
    def __init__(self, EXP, POOL_LENGTH, PCC_POOL):
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL

EDGE=set()
POINT=set()
fnet=open(NETWORK_PATH)
for line in fnet:
    seq=line.rstrip().split('\t')
    if line[0]!='#':
        p1=seq[0]
        p2=seq[1]
        edge=[p1,p2]
        edge.sort()
        edge=':'.join(edge)
        EDGE.add(edge)
        POINT.add(p1)
        POINT.add(p2)
fnet.close()


EXP={}
POOL_LENGTH=0

fpool=open(EXP_POOL_PATH)
for line in fpool:
    try:
        seq=line.rstrip().split('\t')
        tmp=[]
        for one in seq[1:]:
            tmp.append(float(one))
        POOL_LENGTH=len(tmp)
        EXP[seq[0]]=tmp
    except Exception as e:
        pass
fpool.close()

print POOL_LENGTH


from scipy import stats
PCC_POOL={}
for edge in EDGE:
    ps=edge.split(':')
    p1=ps[0]
    p2=ps[1]
    try:
        p1_exp=EXP[p1]
        p2_exp=EXP[p2]
        pcc=stats.pearsonr(p1_exp,p2_exp)[0]
        PCC_POOL[edge]=pcc
    except Exception as e:
        pass


fo=open(OUT_DATA_PATH,'w')
data=Person_PCC_Data( EXP, POOL_LENGTH, PCC_POOL)
pickle.dump(data,fo)
fo.close()




