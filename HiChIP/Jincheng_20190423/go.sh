https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644945
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644946
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644947
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2644948

sudo pip install cooler


#cooler dump -b -t pixels --header --join -r chr3:10M-12M -r2 chr17  GSM2644947_Auxin2days-R1.100000.cool | less
#cooler dump -b -t pixels --join GSM2644945_Untreated-R1.100000.cool > GSM2644945_Untreated-R1.100000.cool.bed

cooler coarsen -k 10 GSM2644945_Untreated-R1.100000.cool -o GSM2644945_Untreated-R1.1000k.cool
cooler dump --join GSM2644945_Untreated-R1.1000k.cool > GSM2644945_Untreated-R1.1000k.cool.bed



if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("HiCcompare")
library(HiCcompare)


hesc100kb <- read.table("GSM2644945_Untreated-R1.100000.cool.bed", header = FALSE)
