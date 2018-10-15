library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

load('Seurat_EXP_mb3076_PCA.Robj')
OUT='exp_LigRec_detail.txt'

RL=read.table('ReceptorLigand.txt.mouse',header=T,sep='\t')

this_data=mb3076@data
this_data=apply(this_data,1,scale)

this_CLUSTER=as.character(mb3076@ident)
CNUM=3


ALLINFO=c('Lig','Rec','LigC','RecC','Score')


rl_row=1
while(rl_row<=length(RL[,1])){
        L=as.character(RL[rl_row,1])
        R=as.character(RL[rl_row,2])
        if(L %in% rownames(EXP@data) & R %in% rownames(EXP@data)){
            Lrow= which(rownames(EXP@data) %in% L)
            Rrow= which(rownames(EXP@data) %in% R)
            Lexp= this_data[Lrow,]
            Rexp= this_data[Rrow,]
            Lscore=c()
            Rscore=c()
            cluster=0
            while(cluster<=(CNUM-1)){
                in_this_cluster=which(this_CLUSTER %in% as.character(cluster))
                lscore=mean(Lexp[in_this_cluster])
                rscore=mean(Rexp[in_this_cluster])
                Lscore=c(Lscore,lscore)
                Rscore=c(Rscore,rscore)
                cluster=cluster+1}
            Lscore[which(is.na(Lscore))]=0
            Rscore[which(is.na(Rscore))]=0
            #Lrank=(rank(Lscore,ties='min')-1)/length(Lscore)
            #Rrank=(rank(Rscore,ties='min')-1)/length(Rscore)
            rr=1
            while(rr<=length(Rscore)){
                cc=1
                while(cc<=length(Lscore)){
                    thisinfo=c(L,R,cc-1,rr-1, (Rscoee[rr]+Lscore[cc])/2)
                    write.table(t(thisinfo),file=OUT,sep='\t',quote=F,row.names=F,col.names=F,append=TRUE)
                    cc=cc+1}
                rr=rr+1}}
        print(rl_row)
        rl_row=rl_row+1}
