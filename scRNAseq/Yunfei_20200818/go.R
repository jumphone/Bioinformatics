
setwd('F:/E_DISK/Project/YunFEI_20200818/DIPG correlation')

source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')


REF_MAT=read.table('TCGA_unified_CORE_ClaNC840.txt',header=T,row.names=1,sep='\t')

REF_MAT=REF_MAT[,2:ncol(REF_MAT)]
TYPE=REF_MAT[1,]
REF_MAT=REF_MAT[2:nrow(REF_MAT),]
RNAME=rownames(REF_MAT)
REF_MAT=apply(REF_MAT,2,as.numeric)
rownames(REF_MAT)=RNAME

TYPE=t(TYPE)
TYPE=cbind(colnames(REF_MAT),TYPE)
LocalRef=.generate_ref(REF_MAT, TYPE, M='mean', min_cell=10)  
#############

SC_MAT_1=read.table('DIPG_KOOLIG2-subtype.txt',header=T,row.names=1,sep='\t')

library(gplots)

out_1=.get_cor(SC_MAT_1, LocalRef, method='spearman',CPU=4, print_step=10)
#out_1=.get_cor(SC_MAT_1, LocalRef, method='pearson',CPU=4, print_step=10)

pdf('ORIG_DIPG_KOOLIG2-subtype.pdf',width=7, height=6)
heatmap.2(out_1,Rowv=T,Colv=T,sepcolor="black", 
          col=colorpanel(256, low="Blue", high="Red", mid="grey85"), 
          scale="col", key=TRUE, keysize=1, symkey=FALSE, 
          densadj = 0.1, density.info="none", trace="none"，
         margins=c(10,12))

dev.off()
sout_1=apply(out_1,2,scale)
rownames(sout_1)=rownames(out_1)
write.table(sout_1,file='ORIG_DIPG_KOOLIG2-subtype.txt',quote=F,row.names=T,col.names=T,sep='\t')

#########
SC_MAT_2=read.table('primary_dipg-rna.tpm-subtype.txt',header=T,row.names=1,sep='\t')

library(gplots)

out_2=.get_cor(SC_MAT_2, LocalRef, method='spearman',CPU=4, print_step=10)
#out_1=.get_cor(SC_MAT_1, LocalRef, method='pearson',CPU=4, print_step=10)

pdf('ORIG_primary_dipg-rna.tpm-subtype.pdf',width=7, height=6)
heatmap.2(out_2,Rowv=T,Colv=T,sepcolor="black", 
          col=colorpanel(256, low="Blue", high="Red", mid="grey85"), 
          scale="col", key=TRUE, keysize=1, symkey=FALSE, 
          densadj = 0.1, density.info="none", trace="none"，
         margins=c(10,12))

dev.off()


sout_2=apply(out_2,2,scale)
rownames(sout_2)=rownames(out_2)
write.table(sout_2,file='ORIG_primary_dipg-rna.tpm-subtype.txt',quote=F,row.names=T,col.names=T,sep='\t')




source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')


SC_MAT=SC_MAT_1
COM_REF_SC=.simple_combine(LocalRef, SC_MAT)$combine
TAG_REF_SC=c(rep('REF',ncol(LocalRef)),rep('SC',ncol(SC_MAT)))
COM_EXP=.combat(COM_REF_SC, TAG_REF_SC)
COM_EXP_REF=COM_EXP[,c(1:ncol(LocalRef))]
COM_EXP_SC=COM_EXP[,c((ncol(LocalRef)+1): ncol(COM_EXP) )]


out_1=.get_cor(COM_EXP_SC, COM_EXP_REF, method='spearman',CPU=4, print_step=10)

pdf('COMBAT_DIPG_KOOLIG2-subtype.pdf',width=7, height=6)
heatmap.2(out_1,Rowv=T,Colv=T,sepcolor="black", 
          col=colorpanel(256, low="Blue", high="Red", mid="grey85"), 
          scale="col", key=TRUE, keysize=1, symkey=FALSE, 
          densadj = 0.1, density.info="none", trace="none"，
         margins=c(10,12))
dev.off()

sout_1=apply(out_1,2,scale)
rownames(sout_1)=rownames(out_1)
write.table(sout_1,file='COMBAT_DIPG_KOOLIG2-subtype.txt',quote=F,row.names=T,col.names=T,sep='\t')



################
SC_MAT=SC_MAT_2
COM_REF_SC=.simple_combine(LocalRef, SC_MAT)$combine
TAG_REF_SC=c(rep('REF',ncol(LocalRef)),rep('SC',ncol(SC_MAT)))
COM_EXP=.combat(COM_REF_SC, TAG_REF_SC)
COM_EXP_REF=COM_EXP[,c(1:ncol(LocalRef))]
COM_EXP_SC=COM_EXP[,c((ncol(LocalRef)+1): ncol(COM_EXP) )]


out_2=.get_cor(COM_EXP_SC, COM_EXP_REF, method='spearman',CPU=4, print_step=10)

pdf('COMBAT_primary_dipg-rna.tpm-subtype.pdf',width=7, height=6)
heatmap.2(out_2,Rowv=T,Colv=T,sepcolor="black", 
          col=colorpanel(256, low="Blue", high="Red", mid="grey85"), 
          scale="col", key=TRUE, keysize=1, symkey=FALSE, 
          densadj = 0.1, density.info="none", trace="none"，
         margins=c(10,12))
dev.off()


sout_2=apply(out_2,2,scale)
rownames(sout_2)=rownames(out_2)
write.table(sout_2,file='COMBAT_primary_dipg-rna.tpm-subtype.txt',quote=F,row.names=T,col.names=T,sep='\t')

























#################################
heatmap.2(x,               #Input必须是matrix
          trace="none",    # trace可以给每个色块中添加一条线，与行平行或者与列平行。其与色块中心的距离代表了这个值被显示的比例。
          scale="none",    # scale在这里
          ColSideColors = plot_color,   # 按照treatment组别给每个subject一个颜色
          dendrogram = "row",   # 生成row的系统发生树
          symbreaks = TRUE,
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)),  # color key, 后面详叙
          breaks = seq(-0.5,0.5,0.05),   # 还是color key
          density.info=c("none"),  # 还是color key
          margins=c(8,16),  # 调整热图大小比例
          cexRow = 0.8, cexCol = 1.0,   # 行列名字体大小
          srtCol = 45, offsetCol = -0.5 # 调整列名的字体倾斜45度，距离热图的距离缩小。
)





