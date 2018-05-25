import sys
fi=open(sys.argv[1])
fo=open('reference_cell.txt','w')
header=fi.readline().rstrip().split('\t')
for one in header:
    if 'NeontalBrain' in one:
        fo.write(one+'\n')
