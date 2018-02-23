fi=open('../data_new/EXP123456.combined.txt')
fo=open('../data_new/EXP123456.combined.txt.cell_color','w')


COL=['begin','red','blue','green','yellow','purple','black']

header=fi.readline().rstrip().split('\t')
lasttag='ok'
i=0
for one in header[1:]:
    tag=one.split('_')[0]
    if tag != lasttag:
        i+=1
        lasttag=tag
    color=COL[i]
    fo.write(one.rstrip()+'\t'+tag+'\t'+color+'\n')
     



