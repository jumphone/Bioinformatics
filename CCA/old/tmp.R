

.dr2adr <- function(DR, B1index, B2index, GROUP, VP, SEED=123){
    set.seed(SEED)
    library(dtw)
    library(MALDIquant)
    library(pcaPP)
    OUT=list()
    OUT$adr=DR
    VALID_PAIR=VP
    ALL_COR=c()   
    ALL_PV=c() 
    ALL_UCOR=c()   
    ALL_UPV=c() 
    index1=B1index
    index2=B2index
  
    vindex1=which(GROUP %in% VP[,1])
    vindex2=which(GROUP %in% VP[,2])
  
    
    print('Start')
    THIS_DR=1
    while(THIS_DR<=ncol(DR)){
        THIS_PC = DR[,THIS_DR]
       
        ########################
        sd_lst1=c()
        mean_lst1=c()
        max_lst1=c()
        min_lst1=c()
        sd_lst2=c()
        mean_lst2=c()  
        max_lst2=c()
        min_lst2=c()
        
        
        
        i=1
        while(i<=nrow(VP)){
            p1=which(GROUP %in% VP[i,1])
            p2=which(GROUP %in% VP[i,2])
            sd1=sd(THIS_PC[p1])
            mean1=mean(THIS_PC[p1])
            sd2=sd(THIS_PC[p2])
            mean2=mean(THIS_PC[p2])
            sd_lst1=c(sd_lst1,sd1)
            sd_lst2=c(sd_lst2,sd2)
            
            max_lst1=c(max_lst1,max(THIS_PC[p1]))
            min_lst1=c(min_lst1,min(THIS_PC[p1]))
            max_lst2=c(max_lst2,max(THIS_PC[p2]))
            min_lst2=c(min_lst2,min(THIS_PC[p2]))
            
            mean_lst1=c(mean_lst1,mean1)
            mean_lst2=c(mean_lst2,mean2)
            i=i+1}
    

        mean_com= apply(cbind(mean_lst1,mean_lst2),1,mean)    
   
      
        .x1_to_com=function(x1){
            #if(x1 <=min_lst1){x1=min_lst1}
            #if(x1 >=max_lst1){x1=max_lst1}
            #if(x1 <=min(mean_lst1)){x1=min(mean_lst1)}
            #if(x1 >=max(mean_lst1)){x1=max(mean_lst1)}
            x1=x1
            dlst1=c()
              
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst1[i]
                this_mean=mean_lst1[i]
                this_d=dnorm(x1,sd=this_sd,mean=this_mean)
                #######################   
                if(x1 > max_lst1[i] | x1 < min_lst1[i]){this_d=0}
                ##################
                if(is.na(this_d)){this_d=0}
                dlst1=c(dlst1,this_d)
                i=i+1} 
            
            if(sum(dlst1)==0){
                    if(x1>max(mean_com){
                    x1=max(mean_com)
                    }else if(x1<min(mean_com){
                    x1=min(mean_com)
                    }else{
                    x1=x1}
                    out=x1
                }else{
                    out=sum(dlst1/sum(dlst1)*mean_com)}
            
            return(out)}
      
        .x2_to_com=function(x2){
            #if(x2 <=min_lst2){x2=min_lst2}
            #if(x2 >=max_lst2){x2=max_lst2}
            #if(x2 <=min(mean_lst2)){x2=min(mean_lst2)}
            #if(x2 >=max(mean_lst2)){x2=max(mean_lst2)}
                
            x2=x2
            dlst2=c()
            
            i=1
            while(i<=nrow(VP)){
                this_sd=sd_lst2[i]
                this_mean=mean_lst2[i]
                this_d=dnorm(x2,sd=this_sd,mean=this_mean)
                #######################   
                if(x2 > max_lst2[i] | x2 < min_lst2[i]){this_d=0}
                ##################
                if(is.na(this_d)){this_d=0}
                dlst2=c(dlst2,this_d)
                i=i+1} 
           
            if(sum(dlst2)==0){
                    if(x2>max(mean_com){
                    x2=max(mean_com)
                    }else if(x2<min(mean_com){
                    x2=min(mean_com)
                    }else{
                    x2=x2}
                    out=x2
                }else{
                    out=sum(dlst2/sum(dlst2)*mean_com)}
            return(out)}
         
         ########################
         #tmp=c(-100:100)
         #tmp1=apply(as.matrix(tmp),1,.x1_to_com)
         #tmp2=apply(as.matrix(tmp),1,.x2_to_com) 
         #plot(tmp1,tmp2)
        
        lst1lst1=apply(as.matrix(DR[index1,THIS_DR]),1,.x1_to_com) 
        lst2lst2=apply(as.matrix(DR[index2,THIS_DR]),1,.x2_to_com)        
       
        OUT$adr[index1,THIS_DR]=lst1lst1
        OUT$adr[index2,THIS_DR]=lst2lst2
        par(mfrow=c(1,2))
        plot(DR[index1,THIS_DR],lst1lst1)
        plot(DR[index2,THIS_DR],lst2lst2)
        
        lst1_mean=c()
        lst2_mean=c()
        i=1
        while(i<=nrow(VALID_PAIR)){
            this_pair=VALID_PAIR[i,]
            this_index1=which(GROUP %in% this_pair[1])
            this_index2=which(GROUP %in% this_pair[2])
            lst1_mean=c(lst1_mean,mean(OUT$adr[this_index1,THIS_DR]))
            lst2_mean=c(lst2_mean,mean(OUT$adr[this_index2,THIS_DR]))
            
            i=i+1}
        
        this_test=cor.test(lst1_mean,lst2_mean)#sum(dist_lst)
        this_ua_test=cor.test(mean_lst1, mean_lst2)
        
        this_cor=this_test$estimate
        this_pv=this_test$p.value
        
        this_un_cor=this_ua_test$estimate
        this_un_pv=this_ua_test$p.value
        
        ALL_COR=c(ALL_COR, this_cor)
        ALL_PV=c(ALL_PV, this_pv) 
        ALL_UCOR=c(ALL_UCOR,this_un_cor)   
        ALL_UPV=c(ALL_UPV,this_un_pv)
        print(THIS_DR)
        
        THIS_DR=THIS_DR+1}
    
    #OUT$cor=ALL_COR
    OUT$cor=ALL_COR
    OUT$pv=ALL_PV
    OUT$ucor=ALL_UCOR
    OUT$upv=ALL_UPV
    print('Finished!!!')
    return(OUT)
   }


DR=bastout$seurat@dr$oldpca@cell.embeddings
GROUP=c(bastout$g1,bastout$g2)
B1index=c(1:length(bastout$g1))
B2index=c((length(bastout$g1)+1):(length(bastout$g1)+length(bastout$g2)))
VP=bastout$vp



OUT=.dr2adr(DR, B1index, B2index, GROUP, VP, SEED=123)
plot(OUT$cor,OUT$ucor)
boxplot(OUT$cor,OUT$ucor)
summary(OUT$cor)
