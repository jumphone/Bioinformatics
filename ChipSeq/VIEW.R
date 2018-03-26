GENE=read.table('../GENE.txt.mouse.EN.bed.sort',sep='\t')
system('mkdir IMAGE')

WT_file='../mapping_mm10/HSC-WT-K27Ac.rmdup.bam.gene.depth'
KO_file='../mapping_mm10/HSC-KO-K27Ac.rmdup.bam.gene.depth'
KI_file='../mapping_mm10/HSC-KI-K27Ac.rmdup.bam.gene.depth'

WT=read.table(WT_file,sep='\t')
KO=read.table(KO_file,sep='\t')
KI=read.table(KI_file,sep='\t')


WTR=2140619/1000000
KOR=8696092/1000000
KIR=5374689/1000000

WT[,3]=WT[,3]/WTR
KO[,3]=KO[,3]/KOR
KI[,3]=KI[,3]/KIR

i=1
while(i<=length(GENE[,1])){
print(i)
chrr=GENE[i,1]
start=GENE[i,2]
end=GENE[i,3]
gene=GENE[i,5]
strand=GENE[i,6]

WT_V=which(WT[,1]==chrr & WT[,2]>=start & WT[,2]<=end)
KO_V=which(KO[,1]==chrr & KO[,2]>=start & KO[,2]<=end)
KI_V=which(KI[,1]==chrr & KI[,2]>=start & KI[,2]<=end)

MAX=0
MAX=max(MAX,WT[WT_V,3])
MAX=max(MAX,KO[KO_V,3])
MAX=max(MAX,KI[KI_V,3])

pdf(paste0('./IMAGE/',gene,'.pdf'),width=20,height=14)
par(mfrow=c(3,1))
plot(WT[WT_V,2],WT[WT_V,3],type='h',xlim=c(start,end),ylim=c(0,MAX),main=paste0(gene,' ',strand),col='blue',xlab='',bty='l',ylab='WT')
plot(KO[KO_V,2],KO[KO_V,3],type='h',xlim=c(start,end),ylim=c(0,MAX),main='',col='red',xlab='',bty='l',ylab='KO')
plot(KI[KI_V,2],KI[KI_V,3],type='h',xlim=c(start,end),ylim=c(0,MAX),main='',col='green',xlab=chrr,bty='l',ylab='KI')

dev.off()
i=i+1}
