
library(circlize)

pdf('CNV.pdf',width=7,height=7)

circos.initializeWithIdeogram(species='mm10')
UPLIMIT=8
LWLIMIT=-2
CEX=0.3
BCOL='grey30'
H=2
##############
bed = read.table('NSC3.BedGraph.bed',sep='\t')
bed[,4][which(bed[,4]>UPLIMIT)]=UPLIMIT
colnames(bed)=c('chr','start','end','value1','value2')

circos.genomicTrackPlotRegion(bed, ylim = c(LWLIMIT, UPLIMIT), panel.fun = function(region, value, ...) {
    #print(value[,2])
   h=H
   cell.xlim = get.cell.meta.data("cell.xlim")
   circos.lines(cell.xlim, c(h, h), col = BCOL)
        
   col = rep('grey50',length(value[,2]))
   col[which(value[,2]==1)]='red'
   col[which(value[,2]==-1)]='blue'
    #circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
   circos.genomicPoints(region, value[,1],col = col, cex = CEX, pch = 16)

    
}, track.height = 0.1)


##############
bed = read.table('NSC4.BedGraph.bed',sep='\t')
bed[,4][which(bed[,4]>UPLIMIT)]=UPLIMIT
colnames(bed)=c('chr','start','end','value1','value2')

circos.genomicTrackPlotRegion(bed, ylim = c(LWLIMIT, UPLIMIT), panel.fun = function(region, value, ...) {
    #print(value[,2])
   h=H
   cell.xlim = get.cell.meta.data("cell.xlim")
   circos.lines(cell.xlim, c(h, h), col = BCOL)
    
   col = rep('grey50',length(value[,2]))
   col[which(value[,2]==1)]='red'
   col[which(value[,2]==-1)]='blue'
    #circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
   circos.genomicPoints(region, value[,1],col = col, cex = CEX, pch = 16)  
}, track.height = 0.1)


##############
bed = read.table('NSC7.BedGraph.bed',sep='\t')
bed[,4][which(bed[,4]>UPLIMIT)]=UPLIMIT
colnames(bed)=c('chr','start','end','value1','value2')

circos.genomicTrackPlotRegion(bed, ylim = c(LWLIMIT, UPLIMIT), panel.fun = function(region, value, ...) {
    #print(value[,2])
   h=H
   cell.xlim = get.cell.meta.data("cell.xlim")
   circos.lines(cell.xlim, c(h, h), col = BCOL)
    
   col = rep('grey50',length(value[,2]))
   col[which(value[,2]==1)]='red'
   col[which(value[,2]==-1)]='blue'
    #circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
   circos.genomicPoints(region, value[,1],col = col, cex = CEX, pch = 16)

    
}, track.height = 0.1)

dev.off()


######################

BW=0.2

########################################################
pdf('CNV_3.pdf',width=7,height=7)

circos.initializeWithIdeogram(species='mm10')
UPLIMIT=8
LWLIMIT=-2
CEX=0.3
BCOL='grey30'
H=2
##############
bed = read.table('NSC3.BedGraph.bed',sep='\t')
bed[,4][which(bed[,4]>UPLIMIT)]=UPLIMIT
colnames(bed)=c('chr','start','end','value1','value2')

circos.genomicTrackPlotRegion(bed, ylim = c(LWLIMIT, UPLIMIT), panel.fun = function(region, value, ...) {
    #print(value[,2])
   h=H
   cell.xlim = get.cell.meta.data("cell.xlim")
   circos.lines(cell.xlim, c(h, h), col = BCOL)
        
   col = rep('grey50',length(value[,2]))
   col[which(value[,2]==1)]='red'
   col[which(value[,2]==-1)]='blue'
    #circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
   circos.genomicPoints(region, value[,1],col = col, cex = CEX, pch = 16)
   
}, track.height = BW)
dev.off()
########################################################



########################################################
pdf('CNV_4.pdf',width=7,height=7)

circos.initializeWithIdeogram(species='mm10')
UPLIMIT=8
LWLIMIT=-2
CEX=0.3
BCOL='grey30'
H=2
##############
bed = read.table('NSC4.BedGraph.bed',sep='\t')
bed[,4][which(bed[,4]>UPLIMIT)]=UPLIMIT
colnames(bed)=c('chr','start','end','value1','value2')

circos.genomicTrackPlotRegion(bed, ylim = c(LWLIMIT, UPLIMIT), panel.fun = function(region, value, ...) {
    #print(value[,2])
   h=H
   cell.xlim = get.cell.meta.data("cell.xlim")
   circos.lines(cell.xlim, c(h, h), col = BCOL)
        
   col = rep('grey50',length(value[,2]))
   col[which(value[,2]==1)]='red'
   col[which(value[,2]==-1)]='blue'
    #circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
   circos.genomicPoints(region, value[,1],col = col, cex = CEX, pch = 16)
   
}, track.height = BW)
dev.off()
########################################################


########################################################
pdf('CNV_7.pdf',width=7,height=7)

circos.initializeWithIdeogram(species='mm10')
UPLIMIT=8
LWLIMIT=-2
CEX=0.3
BCOL='grey30'
H=2
##############
bed = read.table('NSC7.BedGraph.bed',sep='\t')
bed[,4][which(bed[,4]>UPLIMIT)]=UPLIMIT
colnames(bed)=c('chr','start','end','value1','value2')

circos.genomicTrackPlotRegion(bed, ylim = c(LWLIMIT, UPLIMIT), panel.fun = function(region, value, ...) {
    #print(value[,2])
   h=H
   cell.xlim = get.cell.meta.data("cell.xlim")
   circos.lines(cell.xlim, c(h, h), col = BCOL)
        
   col = rep('grey50',length(value[,2]))
   col[which(value[,2]==1)]='red'
   col[which(value[,2]==-1)]='blue'
    #circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
   circos.genomicPoints(region, value[,1],col = col, cex = CEX, pch = 16)
   
}, track.height = BW)
dev.off()
########################################################


