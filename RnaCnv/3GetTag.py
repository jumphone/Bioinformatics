BAND=['red','blue','yellow','green','black','purple','royalblue','indianred']*3000


fi=open('REL_EXP.sorted.txt')
fo=open('cell_tag.txt','w')
line=fi.readline()
seq=line.rstrip().split('\t')

OLD={}
i=0
for one in seq[4:]:
    cell_tag=one.split('_')[0]
    try:
        color=OLD[cell_tag]
    except Exception as e:
        color=BAND[i]
        OLD[cell_tag]=color
        i+=1
    fo.write(cell_tag+'\t'+str(color)+'\n')
OLD={}
i=0
fo=open('gene_tag.txt','w')
for line in fi:
    seq=line.rstrip().split('\t')
    gene_tag=seq[0]
    try:
        color=OLD[gene_tag]
    except Exception as e:
        color=BAND[i]
        OLD[gene_tag]=color
        i+=1
    fo.write(gene_tag+'\t'+str(color)+'\n')
