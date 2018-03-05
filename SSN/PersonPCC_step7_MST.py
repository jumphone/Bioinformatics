import math,sys
import networkx as nx
import matplotlib
matplotlib.use('agg')
import pylab as plt





#ORIGIN=sys.argv[2]#'Mtor'

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
        #if p1==ORIGIN or p2==ORIGIN:
        #    P.add(p1)
        #    P.add(p2)
        w = MAX - score
        G.add_edge(p1,p2,weight=w)
fi.close()

IP=set()
fa = open(sys.argv[1]+'.gene')
for line in fa:
    seq=line.split('\t')
    if int(seq[0])<=500:
        IP.add(seq[1])
    else: 
        break
fa.close()



G_new=nx.Graph()
fi=open(sys.argv[1])
for line in fi:
    seq=line.rstrip().split('\t')
    p1 =seq[0]
    p2=seq[2]
    score=abs(float(seq[4]))
    
    if score > 0: #and ( nx.has_path(G,p1,ORIGIN) or nx.has_path(G,p2,ORIGIN) ):  # ( nx.has_path(G,p1,ORIGIN) or nx.has_path(G,p2,ORIGIN) ) and (p1 in IP and p2 in IP):    #(p1 in P or p2 in P)  #( nx.has_path(G,p1,ORIGIN) or nx.has_path(G,p2,ORIGIN) ):
        w = MAX - score
        G_new.add_edge(p1,p2,weight=w)
fi.close()


MST = nx.minimum_spanning_tree(G_new)
 

#print MST

#fo=open(sys.argv[1]+'.'+ORIGIN+'.MST','w')
fo=open(sys.argv[1]+'.500.MST','w')
fo.write('GENE1\tGENE2\tSCORE\n')
for edge in MST.edges():
    p1=edge[0]
    p2=edge[1]
    if p1 in IP and p2 in IP:
        fo.write(edge[0]+'\t'+edge[1]+'\t1.0\n')

#for one in MST.nodes():
#   print one






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
