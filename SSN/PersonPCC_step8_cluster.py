import sys,subprocess
import networkx as nx

fi=open(sys.argv[1])


G=nx.Graph()
for line in fi:
    seq=line.split('\t')
    if seq[0]!='GENE1':
        G.add_edge(seq[0],seq[1]) 
output=nx.connected_components(G)
ALL=[]
for out in output:
    tmp=[]
    for one in out:
         tmp.append(one)
    ALL.append([len(tmp),tmp])

ALL.sort(reverse=True)
subprocess.Popen("mkdir "+sys.argv[1]+'.cluster',shell=True).wait()
i=0
while i<len(ALL):
    fo=open(sys.argv[1]+'.cluster/cluster.'+str(i+1),'w')
    fo.write('\n'.join(ALL[i][1])+'\n')
    fo.close()
    i=i+1
