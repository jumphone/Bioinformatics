import sys
import cPickle as pickle
STRING_HGNC_PATH = "/home/zhangfeng/disk/project/PersonPCC/database/9606.protein.aliases.v10.5.txt.HGNC_new.symbol"#'/home/zhangfeng/disk/project/PersonPCC/database/9606.protein.aliases.v10.5.txt.HGNC'
NETWORK_PATH = '/home/zhangfeng/disk/project/PersonPCC/database/9606.protein.links.v10.5.txt'
SAMPLE_POOL_PATH = '/home/zhangfeng/disk/project/TCGA/workspace/data_nm_protein_log_FPKM.txt.human.5'



class Person_PCC_Data:
    def __init__(self, string2human , human2string , EXP , POOL_LENGTH , PCC_POOL ):
        self.string2human = string2human
        self.human2string = human2string
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL





string2human={}
human2string={}


fa=open(STRING_HGNC_PATH)
for line in fa:
    seq=line.rstrip().split('\t') 
    string=seq[2]
    human=seq[1]
    string2human[string] = human
    human2string[human] = string
fa.close()


EDGE=set()
POINT=set()
cut_off=700
fnet=open(NETWORK_PATH)
for line in fnet:
    seq=line.rstrip().split(' ')
    try:
        score=float(seq[2])
        if score > cut_off:
            p1=seq[0]
            p2=seq[1]
            edge=[p1,p2]
            edge.sort()
            edge=':'.join(edge)
            EDGE.add(edge)
            POINT.add(p1)
            POINT.add(p2)
    except Exception as e:
        pass
fnet.close()


EXP={}
POOL_LENGTH=0

fpool=open(SAMPLE_POOL_PATH)

for line in fpool:
    seq=line.rstrip().split('\t')
    tmp=[]
    for one in seq[1:]:
        tmp.append(float(one))
    POOL_LENGTH=len(tmp)
    EXP[seq[0]]=tmp

fpool.close()

print POOL_LENGTH


from scipy import stats
PCC_POOL={}
for edge in EDGE:
    ps=edge.split(':')
    p1=ps[0]
    p2=ps[1]
    try:
        p1_human=string2human[p1]
        p2_human=string2human[p2]
        p1_exp=EXP[p1_human]
        p2_exp=EXP[p2_human]
        pcc=stats.pearsonr(p1_exp,p2_exp)[0]
        #if math.isnan(pcc)== False:
        PCC_POOL[edge]=pcc
    except Exception as e:
        pass
import sys
fo=open(sys.argv[1],'w')
data= Person_PCC_Data(string2human , human2string , EXP , POOL_LENGTH , PCC_POOL)
pickle.dump(data,fo)
fo.close()




