setwd('F:/GuoFuKun/NEW1031/GSE114727_RAW/')


.readData=function(PATH,TAG){
    PATH=PATH
    TAG=TAG
    this_data=read.csv(gzfile(PATH), sep=',',header=T)
    this_data=this_data[,c(2:ncol(this_data))]
    this_data=t(this_data)
    this_data[which(is.na(this_data))]=0

    colnames(this_data)=paste0(TAG,'_',c(1:ncol(this_data)))
    print(this_data[1:3,1:3])
    print(dim(this_data))
    return(this_data)

    }


PATHs=c(
        "GSM3148640_BC08_TUMOR3_counts.csv.gz",
        "GSM3148639_BC08_TUMOR2_counts.csv.gz",
        "GSM3148638_BC08_TUMOR1_counts.csv.gz",
        'GSM3148637_BC08_NORMAL1_counts.csv.gz',
        'GSM3148636_BC07_TUMOR4_counts.csv.gz',
        'GSM3148635_BC07_TUMOR3_counts.csv.gz',
        'GSM3148634_BC07_TUMOR2_counts.csv.gz',
        'GSM3148633_BC06_TUMOR3_counts.csv.gz',
        'GSM3148632_BC06_TUMOR2_counts.csv.gz',
        'GSM3148631_BC06_TUMOR1_counts.csv.gz',
        'GSM3148630_BC05_TUMOR4_counts.csv.gz',
        'GSM3148629_BC05_TUMOR3_counts.csv.gz',
        'GSM3148628_BC05_TUMOR2_counts.csv.gz',
        'GSM3148627_BC05_TUMOR1_counts.csv.gz',
        'GSM3148626_BC04_TUMOR6_counts.csv.gz',
        'GSM3148625_BC04_TUMOR5_counts.csv.gz',
        'GSM3148624_BC04_TUMOR4_counts.csv.gz',
        'GSM3148623_BC04_TUMOR3_counts.csv.gz',
        'GSM3148622_BC04_TUMOR2_counts.csv.gz',
        'GSM3148621_BC04_TUMOR1_counts.csv.gz',
        'GSM3148620_BC04_BLOOD7_counts.csv.gz',
        'GSM3148619_BC04_BLOOD6_counts.csv.gz',
        'GSM3148618_BC04_BLOOD5_counts.csv.gz',
        'GSM3148617_BC04_BLOOD4_counts.csv.gz',
        'GSM3148616_BC04_BLOOD3_counts.csv.gz',
        'GSM3148615_BC04_BLOOD2_counts.csv.gz',
        'GSM3148615_BC04_BLOOD2_counts.csv.gz',
        'GSM3148613_BC03_TUMOR5_counts.csv.gz',
        'GSM3148612_BC03_TUMOR3_4_counts.csv.gz',
        'GSM3148611_BC03_TUMOR1_counts.csv.gz',
        'GSM3148610_BC03_NORMAL4_5_counts.csv.gz',
        'GSM3148609_BC03_NORMAL2_3_counts.csv.gz',
        'GSM3148608_BC03_NORMAL1_counts.csv.gz',
        'GSM3148607_BC02_TUMOR4_counts.csv.gz',
        'GSM3148606_BC02_TUMOR3_counts.csv.gz',
        'GSM3148605_BC02_TUMOR2_counts.csv.gz',
        'GSM3148604_BC02_TUMOR1_counts.csv.gz',
        'GSM3148603_BC02_NORMAL3_counts.csv.gz',
        'GSM3148602_BC02_NORMAL2_counts.csv.gz',
        'GSM3148601_BC02_NORMAL1_counts.csv.gz',
        'GSM3148600_BC02_LYMPHNODE6_counts.csv.gz',
        'GSM3148599_BC02_LYMPHNODE5_counts.csv.gz',
        'GSM3148598_BC02_LYMPHNODE4_counts.csv.gz',
        'GSM3148597_BC02_LYMPHNODE3_counts.csv.gz',
        'GSM3148596_BC02_LYMPHNODE2_counts.csv.gz',
        'GSM3148595_BC02_LYMPHNODE1_counts.csv.gz',
        'GSM3148594_BC01_TUMOR4_counts.csv.gz',
        'GSM3148593_BC01_TUMOR3_counts.csv.gz',
        'GSM3148592_BC01_TUMOR2_counts.csv.gz',
        'GSM3148591_BC01_TUMOR1_counts.csv.gz',
        'GSM3148590_BC01_NORMAL4_counts.csv.gz',
        'GSM3148589_BC01_NORMAL3_counts.csv.gz',
        'GSM3148588_BC01_NORMAL2_counts.csv.gz',
        'GSM3148587_BC01_NORMAL1_counts.csv.gz',
        'GSM3148586_BC01_BLOOD3_counts.csv.gz',
        'GSM3148585_BC01_BLOOD1_counts.csv.gz'
        
        )
TAGs=c(  
        "BC08.TUMOR3_counts.csv.gz",
        "BC08.TUMOR2_counts.csv.gz",
        "BC08.TUMOR1_counts.csv.gz",
        'BC08.NORMAL1_counts.csv.gz',
        'BC07.TUMOR4_counts.csv.gz',
        'BC07.TUMOR3_counts.csv.gz',
        'BC07.TUMOR2_counts.csv.gz',
        'BC06.TUMOR3_counts.csv.gz',
        'BC06.TUMOR2_counts.csv.gz',
        'BC06.TUMOR1_counts.csv.gz',
        'BC05.TUMOR4_counts.csv.gz',
        'BC05.TUMOR3_counts.csv.gz',
        'BC05.TUMOR2_counts.csv.gz',
        'BC05.TUMOR1_counts.csv.gz',
        'BC04.TUMOR6_counts.csv.gz',
        'BC04.TUMOR5_counts.csv.gz',
        'BC04.TUMOR4_counts.csv.gz',
        'BC04.TUMOR3_counts.csv.gz',
        'BC04.TUMOR2_counts.csv.gz',
        'BC04.TUMOR1_counts.csv.gz',
        'BC04.BLOOD7_counts.csv.gz',
        'BC04.BLOOD6_counts.csv.gz',
        'BC04.BLOOD5_counts.csv.gz',
        'BC04.BLOOD4_counts.csv.gz',
        'BC04.BLOOD3_counts.csv.gz',
        'BC04.BLOOD2_counts.csv.gz',
        'BC04.BLOOD2_counts.csv.gz',
        'BC03.TUMOR5_counts.csv.gz',
        'BC03.TUMOR3.4_counts.csv.gz',
        'BC03.TUMOR1_counts.csv.gz',
        'BC03.NORMAL4.5_counts.csv.gz',
        'BC03.NORMAL2.3_counts.csv.gz',
        'BC03.NORMAL1_counts.csv.gz',
        'BC02.TUMOR4_counts.csv.gz',
        'BC02.TUMOR3_counts.csv.gz',
        'BC02.TUMOR2_counts.csv.gz',
        'BC02.TUMOR1_counts.csv.gz',
        'BC02.NORMAL3_counts.csv.gz',
        'BC02.NORMAL2_counts.csv.gz',
        'BC02.NORMAL1_counts.csv.gz',
        'BC02.LYMPHNODE6_counts.csv.gz',
        'BC02.LYMPHNODE5_counts.csv.gz',
        'BC02.LYMPHNODE4_counts.csv.gz',
        'BC02.LYMPHNODE3_counts.csv.gz',
        'BC02.LYMPHNODE2_counts.csv.gz',
        'BC02.LYMPHNODE1_counts.csv.gz',
        'BC01.TUMOR4_counts.csv.gz',
        'BC01.TUMOR3_counts.csv.gz',
        'BC01.TUMOR2_counts.csv.gz',
        'BC01.TUMOR1_counts.csv.gz',
        'BC01.NORMAL4_counts.csv.gz',
        'BC01.NORMAL3_counts.csv.gz',
        'BC01.NORMAL2_counts.csv.gz',
        'BC01.NORMAL1_counts.csv.gz',
        'BC01.BLOOD3_counts.csv.gz',
        'BC01.BLOOD1_counts.csv.gz'
     
     )

DATA.LIST=list()

i=1
while(i<=length(PATHs)){
    this_path=PATHs[i]
    this_tag=TAGs[i]
    this_data=.readData(this_path,this_tag)
    DATA.LIST=c(DATA.LIST, list(this_data))
    print(i)
    i=i+1
}

saveRDS(DATA.LIST, 'DATA.LIST.RDS')

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

DATA=DATA.LIST[[1]]
i=2
while(i<=length(DATA.LIST)){
    this_data=DATA.LIST[[i]]
    DATA=.simple_combine(DATA, this_data)$combine

    i=i+1}











