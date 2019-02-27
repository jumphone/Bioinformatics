##########################

from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from MulticoreTSNE import MulticoreTSNE as TSNE
from matplotlib import pyplot as plt
import numpy as np

NT=8
RS=123
PP=30
NC=1
PCP=0.95

##########################

fi=open('GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt')
header=fi.readline()
anno=[]
for line in fi:
    seq=line.rstrip().split('\t')
    anno.append(seq[3])
fi.close()

MSI=list(np.where(np.array(anno)=='MS')[0])
CTI=list(np.where(np.array(anno)=='Ctrl')[0])

##########################
import fileinput

fi=fileinput.input('GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt')
header=fi.readline()
exp=[]
for line in fi:
    exp.append([])
    seq=line.rstrip().split('\t')
    for one in seq[1:]:
        exp[-1].append(float(one))
fi.close()
    
##########################
len(exp) 
exp=np.array(exp)
np.shape(exp)
exp=np.transpose(exp)
np.shape(exp)
#np.save('exp', exp)
#########################
########################

from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from MulticoreTSNE import MulticoreTSNE as TSNE
from matplotlib import pyplot as plt
import numpy as np

NT=8
RS=123
PP=30
NC=1
PCP=0.95
exp=np.load('exp.npy')
###########

fi=open('GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt')
header=fi.readline()
anno=[]
for line in fi:
    seq=line.rstrip().split('\t')
    anno.append(seq[3])
fi.close()

MSI=list(np.where(np.array(anno)=='MS')[0])
CTI=list(np.where(np.array(anno)=='Ctrl')[0])
#################

MSexp=exp[MSI,:]
CTexp=exp[CTI,:]
np.shape(MSexp)
np.shape(CTexp)

#####################
from sklearn import preprocessing

tmp=np.apply_along_axis(np.var, 0, MSexp)
used=np.where(tmp>0)[0]
MSexp=MSexp[:,used]
np.shape(MSexp)

pt = preprocessing.PowerTransformer(method='box-cox', standardize=True)
PT_MSexp=pt.fit_transform((MSexp+1))

PCA_PT_MSexp=PCA(n_components=PCP).fit_transform(PT_MSexp)


MS_embeddings = TSNE(n_jobs=NT,n_components=NC,perplexity=PP,random_state=RS).fit_transform(PCA_PT_MSexp)
MS_x = MS_embeddings[:, 0]

CT_embeddings = TSNE(n_jobs=NT,n_components=NC,perplexity=PP,random_state=RS).fit_transform(CTexp)
CT_x = CT_embeddings[:, 0]

#####################

plt.scatter(vis_x,vis_x,c=digits.target, cmap=plt.cm.get_cmap("jet", 10), marker='.')
plt.colorbar(ticks=range(10))
plt.clim(-0.5, 9.5)
plt.show()


vis_x = embeddings[:, 0]
vis_y = embeddings[:, 1]
plt.scatter(vis_x, vis_y, c=digits.target, cmap=plt.cm.get_cmap("jet", 10), marker='.')
plt.colorbar(ticks=range(10))
plt.clim(-0.5, 9.5)
plt.show()


