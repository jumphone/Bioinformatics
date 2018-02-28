library('plot3D')
phi=0
theta=-45
pdf('3D_PLOT.pdf')
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = phi, theta = theta)
i=0
while(i<=9){
COL=rep('grey90',length(pbmc_stem2_added@pca.rot[,1]))
cc=which(pbmc@ident==i)
COL[cc]='red'
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = phi, theta = theta,col=COL,main=as.character(i))
i=i+1}
dev.off()

