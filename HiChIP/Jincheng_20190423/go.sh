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
con1000kb <- read.table("Untreated.1000k.cool.inter.txt", header = FALSE)
aux1000kb <- read.table("Auxin2days.1000k.cool.inter.txt", header = FALSE)

all.table <- create.hic.table(con1000kb,aux1000kb)


hic.table <- hic_loess(all.table, Plot = TRUE, Plot.smooth = FALSE)

filter_params(hic.table)
hic.table <- hic_compare(hic.table, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)

knitr::kable(head(hic.table))
sig_index=which(hic.table$p.adj<0.05)

sig.hic.table=hic.table[sig_index,]
knitr::kable(head(sig.hic.table))














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



























