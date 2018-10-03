



library(circlize)
a=read.table('SEQ.WT_KO')
CHR=unique(a[,1])
Start=c()
End=c()
for(chr in  CHR){
Start=c(Start,min(a[which(a[,1]==chr),10]))
End=c(End,max(a[which(a[,1]==chr),11]))
}

Example=data.frame(Gene=CHR,Start=as.numeric(Start), End=as.numeric(End))


KO_WT_u_KO = read.table('SEQ.KO_WT_u_KO')
KO_WT_u_WT = read.table('SEQ.KO_WT_u_WT')


circos.genomicInitialize(Example,tickLabelsStartFromZero = TRUE)

bed_list=list(KO_WT_u_KO,KO_WT_u_WT)
KO_WT_col = function(value) {
    col='indianred1'
    if(value[1]==0){col='royalblue1'}

    return(col)
}

circos.genomicTrackPlotRegion(bed_list, ylim=c(0,10), stack = F, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = KO_WT_col(value), border = NA, ...);
    if(value[1]==1){
    circos.text((region[2]-region[1])/2, 5, paste(as.character(value[2]),'%'), facing = "inside", cex = 0.6);
    }

}, bg.border = NA, track.height = 0.1)


