#$1: PersonPCC.data
#$2: Expression
#$3: Output

import cPickle as pickle
import sys
from scipy import stats
print '''
file1: data
file2: EXP, col1: gene, col2: log2(FPKM+1)
file3: Output


'''

class Person_PCC_Data:
    def __init__(self, string2human , human2string , EXP , POOL_LENGTH , PCC_POOL ):
        self.string2human = string2human
        self.human2string = human2string
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL

print 'loading...'
fdata = open(sys.argv[1])
data = pickle.load(fdata)
fdata.close()
print "loading done !"

Person_EXP={}

fi=open(sys.argv[2])
for line in fi:
    seq=line.split('\t')
    Person_EXP[seq[0]]=float(seq[1])
fi.close()

output=[]

fo=open(sys.argv[3],'w')
fo1=open(sys.argv[3]+'.sig','w')
for edge in data.PCC_POOL:
    ps=edge.split(':')
    p1=ps[0]
    p2=ps[1]
    try:
        p1_human=data.string2human[p1]
        p2_human=data.string2human[p2]
        p1_exp_new=data.EXP[p1_human]+[Person_EXP[p1_human]]
        p2_exp_new=data.EXP[p2_human]+[Person_EXP[p2_human]]
        pcc_new = stats.pearsonr(p1_exp_new,p2_exp_new)[0]
        pcc = data.PCC_POOL[edge]
        delta_pcc = pcc_new - pcc            
        z=delta_pcc /( (1-pcc**2)/(data.POOL_LENGTH-1) )
        p=stats.norm.sf(abs(z))*2
        #print p
        line=p1_human+'\t'+p1+'\t'+p2_human+'\t'+p2+'\t'+str(pcc)+'\t'+str(pcc_new)+'\t'+str(z)+'\t'+str(p)
        if '\tnan' not in line and '\t-inf\t' not in line and '\tinf\t' not in line:
            #print p 
            output.append([p,-abs(z),p1,p2,line])
    except Exception as e:
        pass #print e

output.sort()

ALL={}
test_time=len(output)
for one in output:
    pvad=one[0]*test_time
    fo.write(one[-1]+'\t'+str(pvad) +'\n'  )
    if pvad <0.05:
        fo1.write(one[-1]+'\t'+str(pvad) +'\n' )
    p1_= one[2]
    p2 = one[3]
    if p1 in ALL:
            ALL[p1].append(pvad)
    else:
            ALL[p1]=[pvad]
    if p2 in ALL:
            ALL[p2].append(pvad)
    else:
            ALL[p2]=[pvad]

'''
output=[]
for p in ALL:
    num_sigpv=0
    for one in ALL[p]:
        if one < 0.05:
            num_sigpv += 1
    output.append([-num_sigpv,p])
output.sort()
#print output


fo2=open(sys.argv[3]+'.gene','w')
fo3=open(sys.argv[3]+'.gene.sig','w')

for one in output:
    p=one[1]
    human=data.string2human[p]
    #min_pv=str(one[0])
    num_sigpv=str(-one[0])
    
    fo2.write(p+'\t'+human+'\t'+num_sigpv+'\n')
    if float(num_sigpv)>0:
        fo3.write(p+'\t'+human+'\t'+num_sigpv+'\n')

fo.close()
'''
