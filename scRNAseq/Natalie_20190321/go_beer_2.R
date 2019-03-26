
d2=read.table('p100_out_gene_exon_tagged.dge.txt',sep='\t',row.names=1,header=T)


########################

exp_sc_mat=d2
load('Seurat_EXP_cluster.Robj')

pbmc=EXP_cluster
ref_vec=pbmc@dr$tsne@cell.embeddings

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]
ref_tag=cbind(names(pbmc@ident),as.character(pbmc@ident))    
HQRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )   

sc_tag=SCREF(exp_sc_mat, HQRef)$tag2

out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, CPU=4, print_step=10)

########
pdf('OVERLAY.pdf',width=15,height=12)
XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))
plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70')
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red')
dev.off()

saveRDS(out,'overlay.RDS')
saveRDS(sc_tag,'sc_tag.RDS')

