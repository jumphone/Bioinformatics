pseudoExp<-function (pbmc,  gene, ps, clust){
  #ind=which(rownames(pbmc@assays$RNA@data)==gene)
  ind=which(pbmc$cell.type %in% clust)
  ind=rownames(pbmc@meta.data)[ind]
  exprs=pbmc@assays$RNA@data[gene, ind]
  ind2=which(rownames(ps)%in% ind)
  
  ps2=log(ps+1)
  
  all=cbind(ps2[ind2], exprs)
  plot(all[,1],all[,2], xlab="Log Pseudotime",ylab="Normalized Expression")
  c2=pbmc$cell.type[ind]
  for (i in 1:length(c2)){
    if (c2[i]==ct[1]){
      c2[i]="salmon"
    }else if (c2[i]==ct[2]){
      c2[i]="gold2"
    }else if (c2[i]==ct[3]){
      c2[i]="green"
    }else if (c2[i]==ct[4]){
      c2[i]="forestgreen"
    }else if (c2[i]==ct[5]){
      c2[i]="cyan"
    }else if (c2[i]==ct[6]){
      c2[i]="blue"
(all[,1])
  e<-as.matrix(e)
  rownames(e)<-rownames(all)
  ind=which(is.infinite(e))
  e<-e[-ind,1]
  i=cut(e,20)
  g=split(e, i)
  trendMat=matrix(0,20,2)
  for (k in 1:length(g)){
    ids=names(g[[k]])
    xm=mean(as.numeric(all[ids,1]))
    ym=mean(as.numeric(all[ids,2]))
    trendMat[k    }
  }
  all<-cbind(all, c2)
  points(all[,1],all[,2],col=c2, pch=19)
  legend("topright",c(ct[1],ct[2],ct[3],ct[4],ct[5],ct[6]),col=c("salmon","gold2","green","forestgreen","cyan", "blue"),pch=19)
  e=as.numeric,1]=xm
    trendMat[k,2]=ym
  }
  lines(trendMat,col='black', lwd=6)
  title(main=gene)
  
}
