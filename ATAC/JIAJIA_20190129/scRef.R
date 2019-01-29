#library('pcaPP')
########################################
#Author: Feng Zhang
#Email: 15110700005@fudan.edu.cn
#######################################

.get_log_p_sc_given_ref <- function(exp_sc_mat, exp_ref_mat, CPU=4, print_step=10){
    delta = 0.5;
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name   
    ##################
    library(parallel)
    ##################
    Refprob <- function(exp_sc, exp_ref){
    log_p_sc_given_ref = dmultinom(x=exp_sc,log=T,prob=exp_ref)
    return(log_p_sc_given_ref)}
    #################
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    #Step 2. calculate prob
    SINGLE <- function(i){
        exp_sc = as.array(exp_sc_mat[,i])
        log_p_sc_given_ref_list=c()
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            #####
            exp_ref[which(exp_ref==0)]= delta * min(exp_ref[which(exp_ref>0)])
            #####
            log_p_sc_given_ref=Refprob(exp_sc,exp_ref)
            log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
            j=j+1}
        ################################
        if(i%%print_step==1){print(i)}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF = c()
    for(log_p_sc_given_ref_list in RUN){
        LOG_P_SC_GIVEN_REF=cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)}
    #######################################    
    rownames(LOG_P_SC_GIVEN_REF)=colname_ref
    colnames(LOG_P_SC_GIVEN_REF)=colname_sc
    return(LOG_P_SC_GIVEN_REF)
    }

.get_cor  <- function(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10, gene_check=FALSE){
    #method = "pearson", "kendall", "spearman"
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
    ##################
    print('Gene number of exp_sc_mat:')
    print(nrow(exp_sc_mat))
    print('Gene number of exp_ref_mat:')
    print(nrow(exp_ref_mat))
    #################
    library(parallel)
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    ###############
    print('Number of overlapped genes:')
    print(nrow(exp_sc_mat))
    if(gene_check==TRUE){
    print('Press RETURN to continue:')
    scan();}
    ###################
    #Step 2. calculate prob
    SINGLE <- function(i){
        library('pcaPP')
        exp_sc = as.array(exp_sc_mat[,i])
        log_p_sc_given_ref_list=c()
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            #####
            #exp_ref[which(exp_ref==0)]=min(exp_ref[which(exp_ref>0)])
            #####
            #if(method=='rococo'){log_p_sc_given_ref=rococo(exp_sc,exp_ref)} else
            if(method=='kendall'){log_p_sc_given_ref=cor.fk(exp_sc,exp_ref)}
            else{
            log_p_sc_given_ref=cor(exp_sc,exp_ref, method=method)}
            log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
            j=j+1}
        ################################
        if(i%%print_step==1){print(i)}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF = c()
    for(log_p_sc_given_ref_list in RUN){
        LOG_P_SC_GIVEN_REF=cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)}
    #######################################
    rownames(LOG_P_SC_GIVEN_REF)=colname_ref
    colnames(LOG_P_SC_GIVEN_REF)=colname_sc
    return(LOG_P_SC_GIVEN_REF)
    }


.get_p_ref_given_sc <- function(LOG_P_SC_GIVEN_REF){
    #######
    ref2sc <- function(log_p_sc_given_ref_list){
        weight_list = exp(log_p_sc_given_ref_list - max(log_p_sc_given_ref_list))
        p_ref_given_sc = weight_list / sum(weight_list)
        return(p_ref_given_sc)
        }
    #######
    P_REF_GIVEN_SC = apply(LOG_P_SC_GIVEN_REF, 2, ref2sc)
    colnames(P_REF_GIVEN_SC)=colnames(LOG_P_SC_GIVEN_REF)
    rownames(P_REF_GIVEN_SC)=rownames(LOG_P_SC_GIVEN_REF)
    return(P_REF_GIVEN_SC)
    }


.get_tag_max <- function(P_REF_GIVEN_SC){
    RN=rownames(P_REF_GIVEN_SC)
    CN=colnames(P_REF_GIVEN_SC)
    TAG=cbind(CN,rep('NA',length(CN)))
    i=1
    while(i<=length(CN)){
        this_rn_index=which(P_REF_GIVEN_SC[,i] == max(P_REF_GIVEN_SC[,i]))[1]
        TAG[i,2]=RN[this_rn_index]
        i=i+1
        }
    colnames(TAG)=c('cell_id','tag')
    return(TAG)
    }


.get_tag_min <- function(P_REF_GIVEN_SC){
    RN=rownames(P_REF_GIVEN_SC)
    CN=colnames(P_REF_GIVEN_SC)
    TAG=cbind(CN,rep('NA',length(CN)))
    i=1
    while(i<=length(CN)){
        this_rn_index=which(P_REF_GIVEN_SC[,i] == min(P_REF_GIVEN_SC[,i]))[1]
        TAG[i,2]=RN[this_rn_index]
        i=i+1
        }
    colnames(TAG)=c('cell_id','tag')
    return(TAG)
    }


.generate_ref <- function(exp_sc_mat, TAG, min_cell=1, refnames=FALSE){
    NewRef=c()
    TAG[,2]=as.character(TAG[,2])
    if(refnames==FALSE){
        refnames=names(table(TAG[,2]))}
        else{refnames=refnames}
    outnames=c()
    for(one in refnames){
        this_col=which(TAG[,2]==one)
        if(length(this_col)>= min_cell){
            outnames=c(outnames,one)
            if(length(this_col) >1){
                this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
                }
                else{this_new_ref = exp_sc_mat[,this_col]}
            NewRef=cbind(NewRef,this_new_ref)
            }
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }



.compare_two_tag <- function(TAG1, TAG2){
    OUT=c()
    tag1_names=as.character(unique(TAG1[,2]))
    tag2_names=as.character(unique(TAG2[,2]))
    i=1
    while(i<=length(tag1_names)){
        tag1 = tag1_names[i]
        tag1_index = which(TAG1[,2]== tag1)
        j=1
        while(j<=length(tag2_names)){
            tag2 = tag2_names[j] 
            #print(tag2)
            tag2_index = which(TAG2[,2]== tag2)
            over = length(which(tag1_index %in% tag2_index))
            tag1_over = over/length(tag1_index)
            tag2_over = over/length(tag2_index)
            max_over = max(tag1_over, tag2_over)
            OUT=cbind(OUT, c(max_over, tag1, tag1_over ,tag2, tag2_over)) 
            j=j+1
            }
        i=i+1
        }
    OUT=t(OUT)
    #OUT=as.matrix(OUT)
    OUT[,1]=as.numeric(OUT[,1])
    OUT[,3]=as.numeric(OUT[,3])
    OUT[,5]=as.numeric(OUT[,5])
    colnames(OUT)=c('max_over','tag1','tag1_over','tag2','tag2_over')
    OUT=OUT[order(OUT[,1],decreasing=T),]
    return(OUT)
    }




.vec_projection <- function(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec,  method='kendall', nearest_cell=3, random_size=30, random_seed=123, alpha=0.5, min_cell=10, CPU=4, print_step=10){
    delta = 0.5;
    alpha = alpha;
    library(parallel);
    random_seed=random_seed;
    set.seed(random_seed);
    sc_cell_name=colnames(exp_sc_mat);
    tag=sc_tag;
    ref_tag=ref_tag;
    ref_vec=ref_vec;
    exp_ref_mat=exp_ref_mat;
    exp_sc_mat=exp_sc_mat;
    method=method;
    nearest_cell=nearest_cell;
    random_size=random_size;
    min_cell=min_cell;
    CPU=CPU;
    print_step=print_step;
    
    SINGLE = function(i){   
        library('pcaPP')
        Refprob <- function(exp_sc, exp_ref){
        	exp_ref[which(exp_ref==0)] = delta * min(exp_ref[which(exp_ref > 0)]) 
            log_p_sc_given_ref = dmultinom(x=exp_sc,log=T,prob=exp_ref)
            return(log_p_sc_given_ref)}
        .get_dis= function(this_sc, this_ref, method=method){
            exp_sc_mat=this_sc
            exp_ref_mat=this_ref
            exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
            exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
            gene_sc=rownames(exp_sc_mat)
            gene_ref=rownames(exp_ref_mat)
            gene_over= gene_sc[which(gene_sc %in% gene_ref)]
            exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
            exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
            colname_sc=colnames(exp_sc_mat)
            colname_ref=colnames(exp_ref_mat)
            log_p_sc_given_ref_list=c()
            exp_sc = as.array(exp_sc_mat[,1])
            j=1
            while(j<=length(colname_ref)){
                exp_ref = as.array(exp_ref_mat[,j])
                if(method=='multinomial'){
                    log_p_sc_given_ref=Refprob(exp_sc,exp_ref)
                    } else if(method=='kendall'){
                	log_p_sc_given_ref=cor.fk(exp_sc,exp_ref)
                    } else {
                    log_p_sc_given_ref=cor(exp_sc,exp_ref, method=method)
                    }
                
                log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref) 
                j=j+1
                }
            return(log_p_sc_given_ref_list)
            }
        
        this_tag=as.character(tag[i,2])
        #print(this_tag)
        vec_index=which(as.character(ref_tag[,2])==this_tag)
        this_R=min(random_size, length(vec_index))
        vec_index=sample(vec_index, size= this_R, replace = FALSE)         
        this_vec = ref_vec[vec_index,]
        this_ref= exp_ref_mat[,vec_index]
        this_sc = cbind(exp_sc_mat[,i],exp_sc_mat[,i])
        rownames(this_sc) = rownames(exp_sc_mat)
        colnames(this_sc)= c('rep1','rep2')
        this_out = .get_dis(this_sc, this_ref, method=method)
        this_out_rank=rank(-this_out)
        used_index=which(this_out_rank <= nearest_cell)
        
        this_weight = rep(0,length(this_out))
        this_weight[used_index] = 1 #(1-this_out[used_index])/2
        this_weight = this_out_rank * this_weight
        this_weight[used_index] = alpha**this_weight[used_index] 
        this_weight=this_weight/sum(this_weight)

        this_out_vec = t(as.matrix(this_vec)) %*% as.matrix(this_weight)
        this_out_exp='none'
        #this_out_exp = as.matrix(this_ref) %*% as.matrix(this_weight)
        #names(this_out_exp) = rownames(this_ref)
        this_out=list(out_vec=this_out_vec, out_exp=this_out_exp)
        
        if(i%%print_step==1){print(i)}
        return(this_out)
        }
    #windows & unix
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    stopCluster(cl)
    #unix
    #RUN = mclapply(1:length(exp_sc_mat[1,]), SINGLE, mc.cores=CPU)

    OUT_VEC = c()
    #OUT_EXP = c()
    for(this_out in RUN){
        OUT_VEC = cbind(OUT_VEC, this_out$out_vec)
        #OUT_EXP = cbind(OUT_EXP, this_out$out_exp)
        }
    OUT_VEC = t(OUT_VEC)
    rownames(OUT_VEC) = sc_cell_name
    colnames(OUT_VEC) = colnames(ref_vec) 
    #rownames(OUT_EXP) = names(this_out$out_exp)    
    #colnames(OUT_EXP) = sc_cell_name
    #OUT=list(vec=OUT_VEC, exp=OUT_EXP, tag=tag)
    OUT=list(vec=OUT_VEC, tag=tag)
    return(OUT)
    }


SCREF <- function(exp_sc_mat, exp_ref_mat, method1='kendall', method2='multinomial', min_cell=10, CPU=4, print_step=10,gene_check=FALSE){
    print('First-round annotation:')
    print(method1)
    if(method1!='multinomial'){
        out1=.get_cor(exp_sc_mat, exp_ref_mat, method=method1,CPU=CPU, print_step=print_step,gene_check=gene_check)
        } else {
        out1=.get_log_p_sc_given_ref(exp_sc_mat, exp_ref_mat, CPU=CPU, print_step=print_step)
        }
    tag1=.get_tag_max(out1)

    print('Build local reference')
    LocalRef=.generate_ref(exp_sc_mat, tag1, min_cell=min_cell)

    print('Second-round annotation:')
    print(method2)
    if(method2!='multinomial'){
        out2=.get_cor(exp_sc_mat, LocalRef, method=method2,CPU=CPU, print_step=print_step,gene_check=gene_check)
        } else {
        out2=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=CPU, print_step=print_step)
        }
    tag2=.get_tag_max(out2)
    
    output=list()
    output$tag1=tag1
    output$out1=out1
    output$tag2=tag2
    output$out2=out2
    print('Finish!')
    
    return(output)
    }





.trajectory <- function(sim_mat, plot_type='none', random_ratio=0.03, random_seed=123, do.label=TRUE, label_dist=1.2, label_size=3,cell_size=1,plot_size=1.5){

    library(MASS)
    library(ggplot2)
    
    input_value=sim_mat
    random_ratio=random_ratio
    random_seed=random_seed
    label_size=label_size
    plot_size=plot_size
    cell_size=cell_size
    plot_type=plot_type

    set.seed(random_seed)

    CN=length(colnames(input_value))
    N=length(rownames(input_value))

    A=pi/2
    S=2*pi/N
    VEC=c()
    i=1
    while(i<=N){
        A=A+S
        x=cos(A)
        y=sin(A)
        VEC=cbind(VEC,c(x,y))
        i=i+1
    }
    VEC=t(VEC)
    rownames(VEC)=rownames(input_value)
    colnames(VEC)=c('X','Y')

    topN_scale = function(x){
        topN=3
        #s_x=scale(x)
        r=rank(-x,ties.method='random')
        topn=which(r <= topN)
        s_topn=scale(x[topn])
        s_topn[which(s_topn<0)]=0
        y=rep(0,length(x))
        y[topn]=s_topn
        return(y)
        }
    normone = function(x){
        y=x/sum(x)
        return(y)
        }
    
    trascore=function(x){
        y=2*pnorm(x)-1
        y[which(y<0)]=0
        return(y)
        }
    
    tmp = apply(input_value,2,topN_scale)
    #tmp[which(tmp<0)]=0
    tmp = apply(tmp,2,trascore)
    #tmp[which(tmp==0.5)]=0
    tmp = apply(tmp,2,normone)
    colnames(tmp)=colnames(input_value)
    rownames(tmp)=c(rownames(input_value))
    tmp=t(tmp)
    this_vec=tmp %*% VEC 

    r_this_vec=this_vec

    r_this_vec[,1]=this_vec[,1]+rnorm(CN)*random_ratio
    r_this_vec[,2]=this_vec[,2]+rnorm(CN)*random_ratio

    df=data.frame(r_this_vec); colnames(df) = c("x","y")
    
    p=ggplot(data=df,aes(x,y)) +  
      geom_point(colour='grey50',size=cell_size) +
      guides(alpha="none") +
      xlim(-plot_size, plot_size) +
      ylim(-plot_size, plot_size) 
    
    if(plot_type=='polygon'){
        p = p+stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='grey30', contour=T) + 
        scale_fill_continuous(low="green",high="red")  
    }else if(plot_type=='tile'){
        p=p+stat_density2d(aes(alpha=..density..), geom="tile", contour=FALSE)
    }else{p=p}

    seg_vec=c()
    i=1
    while(i<=length(VEC[,1])-1){
    	j=i+1
        while(j <=length(VEC[,1])){
            this_x=VEC[i,1]
            this_y=VEC[i,2]
            this_xe=VEC[j,1]
            this_ye=VEC[j,2]
            seg_vec=cbind(seg_vec, c(this_x,this_y,this_xe,this_ye))
            j=j+1}
        i=i+1}
    seg_vec=t(seg_vec)

    for (i in c(1:length(seg_vec[,1]))) {
        p <- p + geom_segment(x=seg_vec[i,1], y=seg_vec[i,2], xend=seg_vec[i,3], yend=seg_vec[i,4],color="red", linetype="dashed")
        p <- p + geom_point(x=(seg_vec[i,1]+seg_vec[i,3])/2, y=(seg_vec[i,2]+seg_vec[i,4])/2, colour='red', ,size=1.5)
    }

    if(do.label==TRUE){
        target_df=data.frame(VEC* label_dist); colnames(target_df) = c("x","y")
        p=p+geom_text(data=target_df,aes(label=rownames(target_df)), colour='black',size=label_size)
        }

    output=list()
    output$ggplot = p
    output$target_vec=VEC
    output$cell_vec=this_vec
    output$cell_vec_with_random=r_this_vec
    output$mat=tmp
    return(output)
    }



.generate_mst <- function(INPUT, min_cell=10){
    library(igraph)
    INPUT=INPUT
    CNUM=length(INPUT[1,])
    CNAME=colnames(INPUT)

    p1=c()
    p2=c()
    edge_score=c()
    i=1
    while(i<CNUM){
        j=i+1
        while(j<=CNUM){
            p1=c(p1,CNAME[i])
            p2=c(p2,CNAME[j])
            p1_num=length(which(INPUT[,i]>0))
            p2_num=length(which(INPUT[,j]>0))
            p1_and_p2=length(which( INPUT[,i]>0 & INPUT[,j]>0))
            if(p1_and_p2 >= min_cell){
                this_score =  (p1_num * p2_num) / (p1_and_p2)**2}
                else{this_score=NA} 
            edge_score=c(edge_score, this_score)
            j=j+1
        }
        i=i+1
    }

    NO_NA=which(!is.na(edge_score))
    p1=p1[NO_NA]
    p2=p2[NO_NA]
    edge_score=edge_score[NO_NA]
    #beta=10
    #max_tmp=max(edge_score[!is.na(edge_score)])
    #edge_score[is.na(edge_score)] = max_tmp * beta
    NET = cbind(p1,p2) 
    g <- make_graph(t(NET),directed = FALSE)
    MST=mst(g, weights = edge_score, algorithm = NULL)
    return(MST)
    }

#########



.addclust <-function(TAG, TSNE_VEC, MINC=3, MAXC=5, random_seed=123){

    random_seed=random_seed
    TSNE_VEC=TSNE_VEC
    TAG=as.matrix(TAG)
    TAG[,1]=as.character(TAG[,1])
    TAG[,2]=as.character(TAG[,2])
    library(factoextra)
    set.seed(random_seed)
    MAXC=MAXC
    OUT_TAG=TAG
    tag_list=TAG[,2]
    uniq_tag_list=unique(tag_list)
    print('begin')
    i=1
    while(i<=length(uniq_tag_list)){
        this_uniq_tag=uniq_tag_list[i]
        print(this_uniq_tag)
        used_cell=which(tag_list==this_uniq_tag)
        this_tsne_vec=TSNE_VEC[used_cell,]
        this_tag=tag_list[used_cell]

        if(length(this_tag)>=MAXC){
    
            this_opt=fviz_nbclust(this_tsne_vec, kmeans, k.max = MAXC,method = "silhouette")
            this_k = as.numeric(this_opt$data[which(this_opt$data[,2] == max(this_opt$data[,2])),1])
    
            if(this_k<MINC){this_k=1}
            this_out=kmeans(this_tsne_vec,centers=this_k)
    
            this_new_tag=as.character(this_out$cluster)
    
            j=1
            while(j<=length(used_cell)){
                OUT_TAG[used_cell,2][j]=paste0(OUT_TAG[used_cell,2][j],'_',this_new_tag[j])
                j=j+1}
            }
    
        i=i+1;
        }
    print('end')
    return(OUT_TAG)
    }

    
####

.get_gene_cutoff <- function(SYMBOL, exp_sc_mat, BW = FALSE){
    library(pastecs)
    exp_sc_mat=exp_sc_mat
    SYMBOL=SYMBOL
    BW=BW
    GENE=which(rownames(exp_sc_mat)==SYMBOL)
    EXP=as.numeric(exp_sc_mat[GENE,])
    POSEXP=EXP[which(EXP>0)]
    if(BW==FALSE){D=density(POSEXP)}else{D=density(POSEXP,bw=BW)}
    PEAK_PIT=extract(turnpoints(D$y),length(D$y),peak=1,pit=-1)
    LOC=D$x[which(PEAK_PIT==-1)]
    YLIM=c(min(D$y),max(D$y))
    XLIM=c(min(D$x),max(D$x))
    plot(D,main=SYMBOL)
    abline(v=LOC,col='red',lty=3)
    text(x=LOC, y=mean(YLIM), col='red',labels=as.character(round(LOC,2)))
    return(LOC)
} 

.use_gene_cutoff <- function(SYMBOL, exp_sc_mat, cutoff){
    exp_sc_mat=exp_sc_mat
    SYMBOL=SYMBOL
    GENE=which(rownames(exp_sc_mat)==SYMBOL)
    EXP=as.numeric(exp_sc_mat[GENE,])
    TAG=rep(0,length(exp_sc_mat[1,]))
    TAG[which(EXP >=cutoff)]=1
    return(TAG)
}

##2019.01.14#####

.simple_combine <- function(exp_sc_mat1, exp_sc_mat2){    
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
}
    

##2019.01.25####

.check_pos<-function(exp_ref_mat){
    exp_ref_mat=as.matrix(exp_ref_mat)
    exp_ref_mat[which(exp_ref_mat>0)]=1
    SUM=apply(exp_ref_mat,2,sum)
    return(SUM)
}

.trim_pos<-function(exp_ref_mat, trim_number){
    RANK=as.matrix(apply(-exp_ref_mat,2,rank))
    exp_ref_mat=as.matrix(exp_ref_mat)
    exp_ref_mat[which(RANK>trim_number)]=0
    return(exp_ref_mat)
}
    

