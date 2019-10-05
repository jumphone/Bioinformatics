#http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need

#install.packages(c("tm", "SnowballC", "wordcloud", "RColorBrewer", "RCurl", "XML")

                 
#source('https://raw.githubusercontent.com/jumphone/Bioinformatics/master/Wordcloud/rquery_wordcloud.R')
                
#rquery.wordcloud(x=as.character(rnorm(100)), type=c("text"), lang = "english", min.freq = 1,  max.words = 200)


library(wordcloud)

X=c('YES','OK','NO')
F=c(5,4,3)
COL=c('red','blue','black')
set.seed(123)
wordcloud(X,F, min.freq=0, max.words=3,
            random.order=FALSE, rot.per=0.35, 
            use.r.layout=FALSE, colors=COL,ordered.colors=TRUE)



