https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644945
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644946
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644947
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644948

sudo pip install cooler


cooler merge Untreated.cool GSM2644945_Untreated-R1.100000.cool GSM2644946_Untreated-R2.100000.cool
cooler merge Auxin2days.cool GSM2644947_Auxin2days-R1.100000.cool GSM2644948_Auxin2days-R2.100000.cool



cooler coarsen -k 10 Untreated.cool -o Untreated.1000k.cool
cooler coarsen -k 10 Auxin2days.cool -o Auxin2days.1000k.cool


sh cool_inter_mouse.sh Untreated.1000k.cool 
sh cool_inter_mouse.sh Auxin2days.1000k.cool 





#R
#https://bioconductor.org/packages/release/bioc/vignettes/HiCcompare/inst/doc/HiCcompare-vignette.html
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("HiCcompare")


library(HiCcompare)
CHR=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',
'chr16','chr17','chr18','chr19','chrX')
OUT='./OUT/'

for(this_chr in CHR){
#this_chr='chr1'
con1000kb <- read.table(paste0("Untreated.1000k.cool.tmp/",this_chr,'.txt'), header = FALSE)
aux1000kb <- read.table(paste0("Auxin2days.1000k.cool.tmp/",this_chr,'.txt'), header = FALSE)


this.table <- create.hic.table(con1000kb, aux1000kb)
pdf(paste0(OUT,this_chr,'_result.pdf'),width=10,height=7)
hic.table <- hic_loess(this.table, Plot = TRUE, Plot.smooth = FALSE)
filter_params(hic.table)
hic.table <- hic_compare(hic.table, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
dev.off()

#knitr::kable(head(hic.table))
sig_index=which(hic.table$p.value<0.05)

sig.hic.table=hic.table[sig_index,]
#knitr::kable(head(sig.hic.table))
write.table(sig.hic.table,file=paste0(OUT,this_chr,'_sig.txt'),quote=F,sep='\t',row.names=F,col.names=T)
print(this_chr)
}

















