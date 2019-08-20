
setwd('/users/zha8dh/tianlab/HFZ')

D1=read.table(gzfile("./ALL_EXON_DGE/decidua0117_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D2=read.table(gzfile("./ALL_EXON_DGE/decidua0417-2_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D3=read.table(gzfile("./ALL_EXON_DGE/decidua508_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D4=read.table(gzfile("./ALL_EXON_DGE/decidua510_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D5=read.table(gzfile("./ALL_EXON_DGE/decidua514-2_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D6=read.table(gzfile("./ALL_EXON_DGE/decidua-2018_combined_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D7=read.table(gzfile("./ALL_EXON_DGE/decidua-2019_combined_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D8=read.table(gzfile("./ALL_EXON_DGE/decidua20190215_combined_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D9=read.table(gzfile("./ALL_EXON_DGE/Decidua20190420_L3_1000903_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D10=read.table(gzfile("./ALL_EXON_DGE/placenta0423-2_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D11=read.table(gzfile("./ALL_EXON_DGE/placenta508-1_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D12=read.table(gzfile("./ALL_EXON_DGE/placenta508-2_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D13=read.table(gzfile("./ALL_EXON_DGE/placenta514-2_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D14=read.table(gzfile("./ALL_EXON_DGE/placenta-2018_combined_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D15=read.table(gzfile("./ALL_EXON_DGE/placenta-2019_combined_exon_tagged.dge.txt.gz"),header=T,row.names=1)   
D16=read.table(gzfile("./ALL_EXON_DGE/Placenta20190402_L3_1000902_exon_tagged.dge.txt.gz"),header=T,row.names=1)   


colnames(D1)=paste0('decidua0117','_',colnames(D1))
colnames(D2)=paste0('decidua0417.2','_',colnames(D2))
colnames(D3)=paste0('decidua508','_',colnames(D3))
colnames(D4)=paste0('decidua510','_',colnames(D4))
colnames(D5)=paste0('decidua514.2','_',colnames(D5))
colnames(D6)=paste0('decidua2018c','_',colnames(D6))
colnames(D7)=paste0('decidua2019c','_',colnames(D7))
colnames(D8)=paste0('decidua20190215','_',colnames(D8))
colnames(D9)=paste0('decidua20190420','_',colnames(D9))
colnames(D10)=paste0('placenta0423.2','_',colnames(D10))
colnames(D11)=paste0('placenta508.1','_',colnames(D11))
colnames(D12)=paste0('placenta508.2','_',colnames(D12))
colnames(D13)=paste0('placenta514.2','_',colnames(D13))
colnames(D14)=paste0('placenta2018c','_',colnames(D14))
colnames(D15)=paste0('placenta2019c','_',colnames(D15))
colnames(D16)=paste0('placenta20190402','_',colnames(D16))








