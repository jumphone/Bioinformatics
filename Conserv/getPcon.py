import sys,subprocess
fi=open(sys.argv[1])
rangee=int(sys.argv[2])
fo=open(sys.argv[1]+'.Pcon.'+sys.argv[2],'w')
condir='/home/zhangfeng/disk/project/data/annotation/L1/'
#rangee=20
for line in fi:
        seq=line.split('\t')
        chrr=seq[0]
        start=int(seq[1])
        end=int(seq[2])
        inter=[start,end]
        step_tmp=subprocess.Popen(condir+'bigWigToWig -chrom='+chrr+' -start='+str(inter[0])+' -end='+str(inter[1])+' '+'/home/zhangfeng/disk/project/data/annotation/hg19Conservation/hg19.100way.phyloP100way.bw '+condir+'tmpP.wig',shell=True)
        step_tmp.wait()
        ftmp=open(condir+'tmpP.wig')
        fo.write(seq[0]+'\t'+seq[2])
        for line in ftmp:
                if line[0]!='f':
                        fo.write('\t'+line.replace('\n',''))
        fo.write('\n')
        ftmp.close()




fi.close()
fo.close()
