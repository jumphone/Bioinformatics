

load('./MouseSciaticNerve-master/inst/Inj9dBeads/Inj9dBeads.RData')
load('./MouseSciaticNerve-master/inst/Inj9dBeads/Inj9dBeads_savedRes.RData')
plot(data_for_scClustViz$dr_viz,col=as.factor(data_for_scClustViz$cl$res.3.2),pch=16)

plot
TAG=as.numeric(as.character(data_for_scClustViz$cl$res.3.2))
LABEL=rep('NA',length(TAG))

LABEL[which(TAG %in% c(32))]='Mast.Cells'
LABEL[which(TAG %in% c(3,25))]='VSMCs.+.pericytes'
LABEL[which(TAG %in% c(23,31,2,6))]='Monocytes.+.macrophages'
LABEL[which(TAG %in% c(22,18))]='T.cells'
LABEL[which(TAG %in% c(30))]='NK.cells'
LABEL[which(TAG %in% c(9,11,16,19))]='Endothelial.cells'
LABEL[which(TAG %in% c(24))]='Dividing.Schwann.cells'
LABEL[which(TAG %in% c(8,14,15))]='Schwann.cells'
LABEL[which(TAG %in% c(28))]='Dividing.mesenchymal.cells'
LABEL[which(TAG %in% c(13,4,27))]='B.cells'
LABEL[which(TAG %in% c(7,12,26,17,1,5,10,29,20,21))]='Mesenchymal.cells'

source('scRef.R')
this_Ref=.generate_ref(data_for_scClustViz$nge, cbind(LABEL,LABEL), min_cell=1)
saveRDS(this_Ref,file='Inj9dBeads_Ref.RDS')


########

load('./MouseSciaticNerve-master/inst/Inj9dBeadsMesenchymal/Inj9dBeadsMesenchymal.RData')
load('./MouseSciaticNerve-master/inst/Inj9dBeadsMesenchymal/Inj9dBeadsMesenchymal_savedRes.RData')

