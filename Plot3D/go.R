library('plot3D')

STEM2=read.table('STEM3.mouse')
pbmc_stem2=pbmc@data[which(rownames(pbmc@data) %in% STEM2[,3]),]
stem2_var=apply(pbmc_stem2,1,var)
stem2_var_gene=which(rank(-stem2_var,)<=200)

#pbmc_step2_sum=apply(pbmc_stem2[stem2_var_gene,],2,mean)
pbmc_step2_sum=apply(pbmc_stem2[stem2_var_gene,],2,median)
#pbmc_step2_sum=apply(pbmc_stem2,2,quantile,0.75)
pbmc_stem2_added=AddMetaData(object=pbmc,metadata=pbmc_step2_sum,col.name='STEM2')
#VlnPlot(object=pbmc_stem2_added,features.plot=c('STEM2'))

phi=0
theta=60
pdf('3D_PLOT_MEDIAN.pdf')
PCX=1
PCY=2
cor(pbmc_stem2_added@pca.rot[,PCX])
VlnPlot(object=pbmc_stem2_added,features.plot=c('STEM2'))

x=pbmc_stem2_added@pca.rot[,PCX]
y=pbmc_stem2_added@pca.rot[,PCY]
z=pbmc_step2_sum
scatter3D(x,y,z,pch=16,phi = 0, theta = -45)
scatter3D(x,y,z,pch=16,phi = 30, theta = -45)
scatter3D(x,y,z,pch=16,phi = 45, theta = -45)
scatter3D(x,y,z,pch=16,phi = 60, theta = -45)
scatter3D(x,y,z,pch=16,phi = 45, theta = 60)
scatter3D(x,y,z,pch=16,phi = 45, theta = 45)
scatter3D(x,y,z,pch=16,phi = 45, theta = 30)
scatter3D(x,y,z,pch=16,phi = 45, theta = 0)
scatter3D(x,y,z,pch=16,phi = 45, theta = -30)
scatter3D(x,y,z,pch=16,phi = 45, theta = -45)
scatter3D(x,y,z,pch=16,phi = 45, theta = -60)
scatter3D(x,y,z,pch=16,phi = 0, theta = 60)
scatter3D(x,y,z,pch=16,phi = 0, theta = 45)
scatter3D(x,y,z,pch=16,phi = 0, theta = 30)
scatter3D(x,y,z,pch=16,phi = 0, theta = 0)
scatter3D(x,y,z,pch=16,phi = 0, theta = -30)
scatter3D(x,y,z,pch=16,phi = 0, theta = -45)
scatter3D(x,y,z,pch=16,phi = 0, theta = -60)
i=0
while(i<=9){
COL=rep('grey90',length(pbmc_stem2_added@pca.rot[,1]))
cc=which(pbmc@ident==i)
COL[cc]='red'
scatter3D(x=pbmc_stem2_added@pca.rot[,PCX],y=pbmc_stem2_added@pca.rot[,PCY],z=pbmc_step2_sum,pch=16,phi = phi, theta = theta,col=COL,main=as.character(i))
i=i+1}
dev.off()

