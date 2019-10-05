#http://www.sthda.com/english/wiki/word-cloud-generator-in-r-one-killer-function-to-do-everything-you-need

install.packages(c("tm", "SnowballC", "wordcloud", "RColorBrewer", "RCurl", "XML")

rquery.wordcloud(x, type=c("text", "url", "file"), 
        lang="english", excludeWords = NULL, 
        textStemming = FALSE,  colorPalette="Dark2",
        max.words=200)
