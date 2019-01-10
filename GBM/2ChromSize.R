
pc1_pos=read.delim('PC1_POS_OUT.txt',sep='\t',header=T)
pc1_neg=read.delim('PC1_NEG_OUT.txt',sep='\t',header=T)
chromsize=read.table('chromInfo.txt',sep='\t')

CHR=as.character(unique(chromsize[,1]))

pc1_pos_loc=pc1_pos[,c(1:2)]
pc1_neg_loc=pc1_neg[,c(1:2)]

for(this_chr in CHR){
    pc1_pos_loc[which(pc1_pos[,1]==this_chr),2]=pc1_pos_loc[which(pc1_pos[,1]==this_chr),2]/chromsize[which(chromsize[,1]==this_chr),2]
    pc1_neg_loc[which(pc1_neg[,1]==this_chr),2]=pc1_neg_loc[which(pc1_neg[,1]==this_chr),2]/chromsize[which(chromsize[,1]==this_chr),2]
      
}

par(mfrow=c(2,1))
B=50
D=20
X='ChromRelativeLocation'
hist(pc1_pos_loc[,2],breaks=B,density=D, xlab=X,xlim=c(0,1), main='PC1_POS_PEAK')
hist(pc1_neg_loc[,2],breaks=B,density=D, xlab=X,xlim=c(0,1), main='PC1_NEG_PEAK')

ks.test(pc1_pos_loc[,2],pc1_neg_loc[,2])

