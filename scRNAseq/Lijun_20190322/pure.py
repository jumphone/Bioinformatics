fi=open('GSE115011_marton_all_cells.csv')
fo=open('GSE115011_marton_all_cells.csv.pure','w')
old=set()

for line in fi:
    seq=line.rstrip().split(',')
    if seq[0] not in old:
        old.add(seq[0])
        fo.write(line)



