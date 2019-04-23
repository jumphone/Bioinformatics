https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644945
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644946
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644947
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644948

sudo pip install cooler


#cooler dump -b -t pixels --header --join -r chr3:10M-12M -r2 chr17  GSM2644947_Auxin2days-R1.100000.cool | less
#cooler dump -b -t pixels --join GSM2644945_Untreated-R1.100000.cool > GSM2644945_Untreated-R1.100000.cool.bed

cooler coarsen -k 10 GSM2644945_Untreated-R1.100000.cool -o GSM2644945_Untreated-R1.1000k.cool
cooler dump --join -r chr1 GSM2644945_Untreated-R1.1000k.cool > GSM2644945_Untreated-R1.1000k.cool.chr1.txt

cooler coarsen -k 10 GSM2644947_Auxin2days-R1.100000.cool -o GSM2644947_Auxin2days-R1.1000k.cool
cooler dump --join -r chr1  GSM2644947_Auxin2days-R1.1000k.cool > GSM2644947_Auxin2days-R1.1000k.cool.chr1.txt







#R
#https://bioconductor.org/packages/release/bioc/vignettes/HiCcompare/inst/doc/HiCcompare-vignette.html
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("HiCcompare")


library(HiCcompare)
aux1000kb <- read.table("GSM2644947_Auxin2days-R1.1000k.cool.chr1.txt", header = FALSE)
con1000kb <- read.table("GSM2644945_Untreated-R1.1000k.cool.chr1.txt", header = FALSE)

chr1.table <- create.hic.table(aux1000kb, con1000kb, chr = 'chr1')
head(chr1.table)

hic.table <- hic_loess(chr1.table, Plot = TRUE, Plot.smooth = FALSE)

knitr::kable(head(hic.table))
filter_params(hic.table)
hic.table <- hic_compare(hic.table, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)

knitr::kable(head(hic.table))
sig_index=which(hic.table$p.adj<0.05)

sig.hic.table=hic.table[sig_index,]
knitr::kable(head(sig.hic.table))



























