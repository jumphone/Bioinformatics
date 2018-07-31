library(Seurat)

exp_data = read.table('expression.txt',header=T,row.names=1)

EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 0, min.genes=0)

a=read.table('info_changed.txt')

tmp_ident=a[c(3:length(a[,1])),2]
tmp_ident=as.factor(tmp_ident)
names(tmp_ident)=names(EXP@ident)
EXP@ident=tmp_ident
head(EXP@ident)


OPC.PTZ.markers <- FindMarkers(object = EXP, ident.1 = 'OPC.PTZ', thresh.use = 0.25, test.use = "t", only.pos = TRUE)
OPC.markers <- FindMarkers(object = EXP, ident.1 = 'OPC', thresh.use = 0.25, test.use = "t", only.pos = TRUE)
OPC.HigherThan.OPC.PTZ.markers <- FindMarkers(object = EXP, ident.1 = 'OPC',ident.2='OPC.PTZ', thresh.use = 0.25, test.use = "t", only.pos = TRUE)
OPC.PTZ.HigherThan.OPC.markers <- FindMarkers(object = EXP, ident.1 = 'OPC.PTZ',ident.2='OPC', thresh.use = 0.25, test.use = "t", only.pos = TRUE)
write.table(file='OPC.txt',OPC.markers,row.names=T,col.names=T,quote=F,sep='\t')
write.table(file='OPC.PTZ.txt',OPC.PTZ.markers,row.names=T,col.names=T,quote=F,sep='\t')
write.table(file='OPC.HigherThan.OPC.PTZ.txt',OPC.HigherThan.OPC.PTZ.markers,row.names=T,col.names=T,quote=F,sep='\t')
write.table(file='OPC.PTZ.HigherThan.OPC.txt',OPC.PTZ.HigherThan.OPC.markers,row.names=T,col.names=T,quote=F,sep='\t')




GET_MARKER=function(TAG){
    PTZ.markers <- FindMarkers(object = EXP, ident.1 = paste0(TAG,'.PTZ'), thresh.use = 0.25, test.use = "t", only.pos = TRUE)
    NORMAL.markers <- FindMarkers(object = EXP, ident.1 = TAG, thresh.use = 0.25, test.use = "t", only.pos = TRUE)
    NORMAL.HigherThan.PTZ.markers <- FindMarkers(object = EXP, ident.1 = TAG,ident.2=paste0(TAG,'.PTZ'), thresh.use = 0.25, test.use = "t", only.pos = TRUE)
    PTZ.HigherThan.NORMAL.markers <- FindMarkers(object = EXP, ident.1 = paste0(TAG, '.PTZ'),ident.2=TAG, thresh.use = 0.25, test.use = "t", only.pos = TRUE)
    write.table(file=paste0(TAG,'.txt'),NORMAL.markers,row.names=T,col.names=T,quote=F,sep='\t')
    write.table(file=paste0(TAG,'.PTZ.txt'),PTZ.markers,row.names=T,col.names=T,quote=F,sep='\t')
    write.table(file=paste0(TAG,'.HigherThan.',TAG,'.PTZ.txt'), NORMAL.HigherThan.PTZ.markers,row.names=T,col.names=T,quote=F,sep='\t')
    write.table(file=paste0(TAG,'.PTZ','.HigherThan.',TAG,'.txt'),PTZ.HigherThan.NORMAL.markers,row.names=T,col.names=T,quote=F,sep='\t')
}


GET_MARKER('Oligo1')
GET_MARKER('Oligo2')
GET_MARKER('MG')
GET_MARKER('Astro')
