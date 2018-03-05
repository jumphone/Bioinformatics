sudo ln -s /usr/X11 /opt/X11


system("locate libSM.6.dylib")
source("https://bioconductor.org/biocLite.R")
biocLite("metagene")

library(metagene)
library(ggplot2)

bam_files=c('H3K27Ac_OPCs.rmdup.bam')
regions=c('H3K27Ac_OPCs_peaks_TypicalEnhancers.bed')
mg <- metagene$new(regions=regions,bam_file=bam_files,padding_size=1500)
mg$produce_data_frame()
DF=mg$get_data_frame()

pdf('demo.pdf',width=20,height=10)
plot_metagene(DF) +  scale_x_continuous(limits=c(0,100)) +scale_y_continuous(limits=c(0,5))
plot_metagene(DF) 
dev.off()
