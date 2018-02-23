library(methylKit)
setwd('XXX/Methy/data/SAM/')
#setwd('XXX/Methy/data/REFSAM/')
QUAL=0
OUT='XXX/Methy/data/METHY'
file_name='_1_val_1_bismark_bt2_pe.sam.sort.sam'
TYPE='CpG'
REF='hg19'


TAG='AO2'
processBismarkAln(location=c(paste0(TAG,file_name)), sample.id=c(TAG),assembly=REF,read.context=TYPE,save.folder=OUT,mincov = 0,minqual = QUAL, phred64 = FALSE,)
TAG='DIPG13'
processBismarkAln(location=c(paste0(TAG,file_name)), sample.id=c(TAG),assembly=REF,read.context=TYPE,save.folder=OUT,mincov = 0,minqual = QUAL, phred64 = FALSE,)
TAG='DIPG4'
processBismarkAln(location=c(paste0(TAG,file_name)), sample.id=c(TAG),assembly=REF,read.context=TYPE,save.folder=OUT,mincov = 0,minqual = QUAL, phred64 = FALSE,)
TAG='DIPG6'
processBismarkAln(location=c(paste0(TAG,file_name)), sample.id=c(TAG),assembly=REF,read.context=TYPE,save.folder=OUT,mincov = 0,minqual = QUAL, phred64 = FALSE,)
TAG='DIPGC1'
processBismarkAln(location=c(paste0(TAG,file_name)), sample.id=c(TAG),assembly=REF,read.context=TYPE,save.folder=OUT,mincov = 0,minqual = QUAL, phred64 = FALSE,)
TAG='DIPGC2'
processBismarkAln(location=c(paste0(TAG,file_name)), sample.id=c(TAG),assembly=REF,read.context=TYPE,save.folder=OUT,mincov = 0,minqual = QUAL, phred64 = FALSE,)


setwd('XXX/Methy/data/METHY/')

file.list=list("AO2_CpG.txt","DIPG13_CpG.txt","DIPG4_CpG.txt", "DIPG6_CpG.txt","DIPGC1_CpG.txt","DIPGC2_CpG.txt")
myobj=methRead(file.list,
           sample.id=list("AO2","DIPG13","DIPG4","DIPG6","DIPGC1","DIPGC2"),
           treatment=c(0,0,0,0,1,0),
           assembly="hg19",
           pipeline='bismark',
           mincov=1,
           context="CpG"
           )


pdf('PLOT.pdf')
i=1
while(i <=6){
getMethylationStats(myobj[[i]],plot=T,both.strands=FALSE)
getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
i=i+1
}
dev.off()

pdf('RESULT.pdf')
filtered.myobj=filterByCoverage(myobj,lo.count=NULL,lo.perc=NULL,hi.count=NULL,hi.perc=NULL)
meth=unite(filtered.myobj, destrand=FALSE)
getCorrelation(meth,plot=TRUE)
#getCorrelation(meth,plot=TRUE,method='spearman')
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)
dev.off()


      
      
      
           
           

