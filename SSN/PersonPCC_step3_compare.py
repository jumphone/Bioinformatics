import sys
#compared z-value =  $2 - $1

print '''
file1: step1 output 1
file2: step1 output 2
output: 2-1


'''


PCC_origin={}
PCC1={}
PCC2={}

PAIR1={}
f1=open(sys.argv[1])
 
for line in f1:
    seq=line.rstrip().split('\t')
    tag=':'.join(seq[0:4])
    #print seq
#    try:
    z=float(seq[6])
    PCC1[tag]=seq[5]
    PCC_origin[tag]=seq[4]
 #   except Exception as e:
  #      print seq
    PAIR1[tag]=z


PAIR2={}
    
f2=open(sys.argv[2])

for line in f2:
    seq=line.rstrip().split('\t')
    tag=':'.join(seq[0:4])
    z=float(seq[6])
    PAIR2[tag]=z
    PCC2[tag]=seq[5]



output=[]
for one in PAIR1:
    if one in PAIR2:
        z=PAIR2[one]-PAIR1[one]
        pc0=PCC_origin[one]
        pc1=PCC1[one]
        pc2=PCC2[one]
        
        tag=one.replace(':', '\t')
        output.append([abs(z),z,tag,pc0,pc1,pc2])


output.sort(reverse=True) 

fo=open(sys.argv[3],'w')

for one in output:
    #if one[0] > 1.97:
        fo.write(one[2]+'\t'+str(one[1])+'\t'+str(one[3])+'\t'+str(one[4])+'\t'+str(one[5])+'\n')


    
