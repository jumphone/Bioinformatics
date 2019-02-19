library('Seurat')
source('scRef.R')

####################
load('Feng1922WorkspaceSubmission.RData')
pbmc=mmap
DimPlot(pbmc, reduction.use='umap')
 
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]


####################
exp_ref_mat=read.table('GSE123025_exp.txt',header=T,row.names=1,sep='\t',check.name=F)
TAG=read.table('GSE123025_label.txt',header=T,sep='\t',check.name=F)
MAP=read.table('MAP.txt',header=T,sep='\t')


MAP[,1]=as.character(MAP[,1])
MAP[,2]=as.character(MAP[,2])
TAG[,15]=as.character(TAG[,15])

tmp=c()
i=1
while(i<=length(TAG[,15])){
   if(TAG[i,15] %in% MAP[,1]){
       this_tmp=MAP[which(MAP[,1]==TAG[i,15]),2]
   }else{this_tmp='NA'}
   tmp=c(tmp,this_tmp)
   i=i+1
}

ref_tag=cbind(as.character(TAG[,1]),tmp)

V=which(tmp!='NA')

exp_ref_mat=exp_ref_mat[,V]
ref_tag=ref_tag[V,]
####################


####################
pbmc=CreateSeuratObject(raw.data = exp_ref_mat, min.cells = 0, min.genes = 0, project = "10X_PBMC")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
pbmc <- RunPCA(object = pbmc, pcs.compute=50,pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
PCElbowPlot(object = pbmc,num.pc=50)
pbmc <- RunUMAP(pbmc, dims.use=1:20)
pbmc@meta.data$tag=ref_tag[,2]
DimPlot(pbmc, reduction.use = "umap",pt.size=0.5,group.by='tag')
######################

ref_vec=pbmc@dr$umap@cell.embeddings
#######################

LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )   
sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2

out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=1000, 
        random_seed=123, CPU=4, print_step=10)



mmap@meta.data$tag=sc_tag[,2]

pdf('ALL.pdf',width=12,height=10)

DimPlot(pbmc, reduction.use = "umap",pt.size=0.5,group.by='tag')
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70')
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red')
pie(table(sc_tag[,2]))
DimPlot(mmap, reduction.use = "umap",do.label=T,pt.size=0.5,group.by='tag')
dev.off()



save.image(file='ALL.RData')

