library('plot3D')

STEM2=read.table('STEM2.mouse')
pbmc_stem2=pbmc@data[which(rownames(pbmc@data) %in% STEM2[,3]),]
stem2_var=apply(pbmc_stem2,1,var)
stem2_var_gene=which(rank(-stem2_var,)<=10)

pbmc_step2_sum=apply(pbmc_stem2[stem2_var_gene,],2,mean)
#pbmc_step2_sum=apply(pbmc_stem2,2,quantile,0.9)
pbmc_stem2_added=AddMetaData(object=pbmc,metadata=pbmc_step2_sum,col.name='STEM2')
VlnPlot(object=pbmc_stem2_added,features.plot=c('STEM2'))

phi=0
theta=-45
pdf('3D_PLOT.pdf')
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = -45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 30, theta = -45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = -45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 60, theta = -45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = 60)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = 45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = 30)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = 0)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = -30)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = -45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 45, theta = -60)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = 60)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = 45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = 30)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = 0)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = -30)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = -45)
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = 0, theta = -60)
i=0
while(i<=9){
COL=rep('grey90',length(pbmc_stem2_added@pca.rot[,1]))
cc=which(pbmc@ident==i)
COL[cc]='red'
scatter3D(x=pbmc_stem2_added@pca.rot[,1],y=pbmc_stem2_added@pca.rot[,2],z=pbmc_step2_sum,pch=16,phi = phi, theta = theta,col=COL,main=as.character(i))
i=i+1}
dev.off()

