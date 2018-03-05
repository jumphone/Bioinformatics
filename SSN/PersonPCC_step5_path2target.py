import math,sys
import networkx as nx
import matplotlib
matplotlib.use('agg')
import pylab as plt





ORIGIN=sys.argv[2]#'Mtor'

MAX=100000



fi=open(sys.argv[1])
G=nx.Graph()
P=set()
for line in fi:
    seq=line.rstrip().split('\t') 
    p1=seq[0]
    p2=seq[2]
    score=abs(float(seq[4]))
    if score >0:
        P.add(p1)
        P.add(p2)
        w =  MAX - score             #math.log(abs(float(seq[4]))+1,2)
        G.add_edge(p1,p2,weight=w)
fi.close()

OP=set()
fi=open(sys.argv[1]+'.gene')

limit=int(sys.argv[3])

fo=open(sys.argv[1]+'.gene.path2'+ORIGIN+'.'+str(limit),'w')

#i=1
i=0
for line in fi:
    if i >=limit:
        break 
    seq=line.rstrip().split('\t')
    p1=ORIGIN
    p2=seq[1]
#    print str(i)+'\r',;i=i+1   
    path='None'
    try: 
        path=','.join(nx.shortest_path(G,p1,p2,weight="weight"))
    except Exception as e:
        pass
    fo.write(p1+'\t'+p2+'\t'+path+'\n')
    i=i+1
'''
fi=open(sys.argv[1])
G_new=nx.Graph()
for line in fi:
    seq=line.rstrip().split('\t')
    p1=seq[0]
    p2=seq[2]
    if p1 in OP and p2 in OP:
        score=abs(float(seq[4]))
        w =  1 / (score+1)#math.log(abs(float(seq[4]))+1,2)
        G_new.add_edge(p1,p2,weight=w)

fi.close()
print "loaded !!!"
#MST=nx.minimum_spanning_tree(G_new)
'''







'''

#pos=nx.spring_layout(G_new)
pos=nx.circular_layout(G_new)
pos_new={}
import numpy
for one in pos:
    pos_new[one]= numpy.asarray([pos[one][0]*10000,pos[one][1]*10000])
#print pos
fig = plt.figure()
fig.figsize=(10000,10000)
dpi=1000
ax=fig.add_subplot(1, 1, 1, frameon=False)
ax.set_xlim(-15000, 15000)
ax.set_ylim(-15000, 15000)
#nx.draw_networkx(MST,with_label=True)
nx.draw_networkx(G_new,pos_new,with_label=True,ax=ax)
plt.savefig('ok.png')
'''
#print MST.nodes(data=True)
#print G.nodes(data=True)
