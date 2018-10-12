gene_exp <-as.matrix(read.table(file=“abc.txt",header=T, row.names=1))
gene_exp_s=t(apply(gene_exp,1,scale,scale=F))
rownames(gene_exp_s)=rownames(gene_exp)
colnames(gene_exp_s)=colnames(gene_exp)
RR=gene_exp_s[,2]-gene_exp_s[,1]
gene_exp_s=gene_exp_s[order(RR),] 
gene_exp_s[which(gene_exp_s[,1]>15),1]=15 
gene_exp_s[which(gene_exp_s[,2]>15),2]=15
gene_exp_s[which(gene_exp_s[,1]< -15),1]= -15
gene_exp_s[which(gene_exp_s[,2]< -15),2]= -15
pdf(“abc.pdf",6,10)
heatmap.2(gene_exp_s,Rowv=NA,Colv=NA,sepcolor="black", col=colorpanel(256, low="Blue", high="Red", mid="White"), scale="none", key=TRUE, keysize=1, symkey=FALSE, densadj = 0.1, density.info="none", trace="none")
dev.off()
