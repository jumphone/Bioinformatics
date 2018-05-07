import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.tsv','w')
fo.write(fi.read().replace('\r','\n'))
