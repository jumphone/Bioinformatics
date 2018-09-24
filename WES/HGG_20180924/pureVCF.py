import sys


def convertinfo(loc,  ref, alt):
    loc=int(loc)
    ref=ref.upper()
    alt=alt.upper()
    if len(ref)==1 and len(alt)==1:
        return([[loc,ref,alt]])
    elif len(ref)>1 and len(alt)==1:
        return([[loc,ref,alt]])
    elif len(ref)==1 and len(alt)>1:
        return([[loc,ref,alt]])
    elif len(ref)>1 and len(alt)>1:
        out=[]
        if len(ref)==len(alt):
            i=0
            while i<len(ref):
                if ref[i]!=alt[i]:
                   out.append([loc+i,ref[i],alt[i]])
                i=i+1
        elif len(ref) > len(alt):
            i=1
            while i<len(alt):
                if alt[:i] == ref[:i] and alt[ (len(alt)-len(alt)+i):] == ref[(len(ref)-len(alt)+i):]:
                    out.append( [loc+i-1, ref[(i-1):(len(ref)-len(alt)+i)], alt[i-1]] )
                    i=len(alt)
                i=i+1
        elif len(alt) > len(ref):
            i=1
            while i<len(ref):
                if alt[:i] == ref[:i] and ref[ (len(ref)-len(ref)+i):] == alt[(len(alt)-len(ref)+i):]:
                    out.append( [loc+i-1, ref[i-1],alt[(i-1):(len(alt)-len(ref)+i)]] )
                    i=len(ref)
                i=i+1
        return out


fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.pure','w')
for line in fi:
    if line[0]=='#':
        fo.write(line)
    else:

        seq=line.rstrip().split('\t')
        CHR=seq[0]
        LOC=int(seq[1])
        REF=seq[3]
        ALT=seq[4]
        ALT=ALT.replace(',<*>','')
        INFO=seq[5:]

        MR=0
        try:
            TTT=seq[8].split(':')
            III=seq[9].split(':')
            AD=float(III[TTT.index('AD')].split(',')[1])
            DP=float(III[TTT.index('DP')].split(',')[0])
            MR=AD/DP
        except Exception as e:
            #print line
            pass

        if MR >=0.1:
            OK=convertinfo(LOC,REF,ALT)
            if ',' not in ALT and len(OK)>0:

                for one in OK:
                    out=[CHR,str(one[0]),'.',one[1],one[2]]+INFO
                    fo.write('\t'.join(out)+'\n')
            else:
                fo.write(line)


fo.close()
