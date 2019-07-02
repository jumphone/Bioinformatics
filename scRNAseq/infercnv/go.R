library(Seurat)
library(infercnv)

###############################################
this_study='BT309'
pbmc=readRDS('BT309_pbmc_tutorial.rds')
TYPE=as.character(Idents(pbmc))
REF_TYPE=c('Microglia','Pericyte','Photoreceptor cells')

##############################################
names(TYPE)=colnames(pbmc)
AL=TYPE

print('Start')
print(this_study)    
this_type=AL
this_exp=as.matrix(pbmc@assays$RNA@data)
    
this_anno=cbind(colnames(this_exp),this_type)
this_anno_file=paste0(this_study,'.anno.txt')
write.table(this_anno,file=this_anno_file,sep='\t',row.names=F,col.names=F,quote=F)


this_cnv = CreateInfercnvObject(raw_counts_matrix=this_exp,
                                    annotations_file= this_anno_file,
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c(REF_TYPE) )
                                    
this_cnv = infercnv::run(this_cnv,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0('./CNV/',this_study),  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             )

    
saveRDS(this_type, paste0('./CNV/',this_study,'/this_type.RDS'))
saveRDS(this_exp, paste0('./CNV/',this_study,'/this_exp.RDS'))
saveRDS(this_cnv, paste0('./CNV/',this_study,'/this_cnv.RDS'))




