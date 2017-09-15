library(tm)
library(SnowballC)
library(slam)
library(readr)
library(data.table)
library(dplyr)
# con = file("sub_twitter.txt", encoding = 'UTF-8')
# sub_twitter <- read_lines(con)
# con = file("sub_news.txt", encoding = 'UTF-8')
# sub_news <- read_lines(con)
# con = file("sub_blogs.txt", encoding = 'UTF-8')
# sub_blogs <- read_lines(con)
con = file("sub_sample.txt", encoding = 'UTF-8')
sample <- read_lines(con)
addTips = function(doc) {
    paste('<s>', doc, '</s>')
}
CrtCleanVCorpus<-function(vector)
{corpus<-VCorpus(VectorSource(vector), 
                 readerControl = list(language = "lat",
                                      load = TRUE)    )
corpus <- tm_map(corpus, removePunctuation)
corpus<-tm_map(corpus, content_transformer(tolower))
corpus<-tm_map(corpus, removeNumbers)
#corpus<-tm_map(corpus, stemDocument)  
corpus<-tm_map(corpus, stripWhitespace)
corpus <- tm_map(corpus, content_transformer(addTips))
}

#generate Corpus objective
data<-CrtCleanVCorpus(sample[1:250000]) 
testData<-CrtCleanVCorpus(sample[250001:250100]) #size of test sample maybe changed
#generate uni-gram 
dtm = DocumentTermMatrix(data, control=list(wordLengths=c(1,Inf)))
dtm = dtm[slam::row_sums(dtm)>0,] #remove empty documents
#control for constructing 2-grams tdm 
BigramTokenizer<- function(x)
        {unlist(lapply(ngrams(words(x), 2), paste, collapse = " "), 
         use.names = FALSE) }
ctrl = list (tokenize = BigramTokenizer, wordLengths=c(1,Inf))
bi_dtm = DocumentTermMatrix(data, control = ctrl)
bi_dtm = bi_dtm[slam::row_sums(bi_dtm)>0,] #remove documents with zero bigrams
#control for constructing 3-grams tdm 
TrigramTokenizer<-function(x)
        {unlist(lapply(ngrams(words(x), 3), paste, collapse = " "), 
         use.names = FALSE) }
ctrl = list (tokenize = TrigramTokenizer)
tri_dtm <- DocumentTermMatrix(data, control = ctrl)
tri_dtm = tri_dtm[slam::row_sums(tri_dtm)>0,] 
#control for constructing 4-grams tdm 
four_gramTokenizer<-function(x)
{unlist(lapply(ngrams(words(x), 4), paste, collapse = " "), 
        use.names = FALSE) }
ctrl = list (tokenize = four_gramTokenizer)
four_dtm <- DocumentTermMatrix(data, control = ctrl)
four_dtm = four_dtm[slam::row_sums(four_dtm)>0,]
#function to generate term-frequncy table 
freqTable<-function(dtm, N = ncol(dtm))
{
    m = slam::col_sums(dtm, na.rm = T) #col sum up the dtm to get the term counts
    v = sort(m, decreasing=TRUE)
    df = data.table(term = as.character(names(v)), freq = v)
    return(df[1:N,])
}

termFreq<-freqTable(dtm = dtm)
termFreq$term = as.character(termFreq$term)
bi_termFreq<-freqTable(dtm = bi_dtm)
bi_termFreq$term = as.character(bi_termFreq$term)
tri_termFreq<-freqTable(dtm = tri_dtm)
tri_termFreq$term = as.character(tri_termFreq$term)
four_termFreq<-freqTable(dtm = four_dtm)
four_termFreq$term = as.character(four_termFreq$term)

write.csv(termFreq, 'uni_gram_freq.csv', row.names = F)
write.csv(bi_termFreq, 'bi_gram_freq.csv', row.names = F)
write.csv(tri_termFreq, 'tri_gram_freq.csv', row.names = F)
write.csv(four_termFreq, 'four_gram_freq.csv', row.names = F)

df1 = fread("uni_gram_freq.csv")
df2 = fread("bi_gram_freq.csv")
df3 = fread("tri_gram_freq.csv")
df4 = fread("four_gram_freq.csv")
###############################################################################################################
#df3 is a data.table with colomns terms and freq. Now the function getFirstTerms/
#getLastTerms split df3 to threeGramTable containing FirstTerm, lastTerm, freq
#we do the same thing for 2-gram and 4gram table.
getFirstTerms = function(char, n) {
    wordList = unlist(strsplit(char, " ", fixed = T))
    FirstTerms = paste(c(wordList[1:n]), collapse = ' ')
    return(FirstTerms)
}
getLastTerms = function(char, n) {
    wordList = unlist(strsplit(char, " ", fixed = T))
    LastTerms=paste(c(wordList[(length(wordList)-n+1):length(wordList)]),collapse = ' ')
    return(LastTerms)
}
oneGramTable = df1
colnames(oneGramTable) = c('firstTerms', 'freq')
oneGramTable=oneGramTable[order(firstTerms, -freq) ]

threeGramTable = df3
threeGramTable$firstTerms = unlist(lapply(df3$term, getFirstTerms, n = 2))
threeGramTable$lastTerm = unlist(lapply(df3$term, getLastTerms, n = 1))
threeGramTable=threeGramTable[with(threeGramTable, order(firstTerms, -freq)), ]
threeGramTable <- threeGramTable[, list(firstTerms, lastTerm, freq)]

twoGramTable = df2
twoGramTable$firstTerms = unlist(lapply(df2$term, getFirstTerms, n = 1))
twoGramTable$lastTerm = unlist(lapply(df2$term, getLastTerms, n = 1))
twoGramTable=twoGramTable[with(twoGramTable, order(firstTerms, -freq)), ]
twoGramTable <- twoGramTable[, list(firstTerms, lastTerm, freq)]

fourGramTable = df4
fourGramTable$firstTerms = unlist(lapply(df4$term, getFirstTerms, n = 3))
fourGramTable$lastTerm = unlist(lapply(df4$term, getLastTerms, n = 1))
fourGramTable=fourGramTable[with(fourGramTable, order(firstTerms, -freq)), ]
fourGramTable <- fourGramTable[, list(firstTerms, lastTerm, freq)]

#remove the wired rows, should be fixed later
oneGramTable = oneGramTable[6:nrow(oneGramTable)]
twoGramTable = twoGramTable[6:nrow(twoGramTable)]
threeGramTable = threeGramTable[4:nrow(threeGramTable)]
fourGramTable = fourGramTable[4:nrow(fourGramTable)]

write.csv(oneGramTable, 'oneGramTable.csv', row.names = F)
write.csv(twoGramTable, 'twoGramTable.csv', row.names = F)
write.csv(threeGramTable, 'threeGramTable.csv', row.names = F)
write.csv(fourGramTable, 'fourGramTable.csv', row.names = F)


#start here when model is fixed
oneGramTable = fread("oneGramTable.csv")
twoGramTable = fread("twoGramTable.csv")
threeGramTable = fread("threeGramTable.csv")
fourGramTable = fread("fourGramTable.csv")

##########################################################################################################
# Now we make extended 3-gram table with Discount column and Remaining Probabilities.
# Apply the formula:
# d = r* / r = ((r+1)/r)(n_(r+1)/n_r)
# For example, for frequency = 5, d = (6/5) * (N6/N5)
# leftOverProb = Beta, i.e. contribution of a m-gram (count of C(W_1:m))
# to the probabiliy of un-seem m-grams

calcLeftOverProb = function(freq, discount){
    all_freq = sum(freq)
    
    return(1-sum((discount*freq)/all_freq))
}
#add a column good-turning discount coef. to ngramTable
goodTurningDiscountCoef = function (nGramTable, k) {
    # Supposed table "nGramTable" as above, we want to add a "discount" column
    # if k>0 for four-gram and >5 for three-gram discount coef. = 1.0
    nGramTable$discount = rep(1, nrow(nGramTable))
    for(i in k:1) {
        currRTimes = i
        nextRTimes = currRTimes + 1
        
        currN = nrow(nGramTable[freq == currRTimes])
        nextN = nrow(nGramTable[freq == nextRTimes])
        
        currd = (nextRTimes / currRTimes) * (nextN / currN) # assumption: 0 < d < 1
        
        #specific use of data.table: for rows with this condition, add column discount as such
        nGramTable[freq == currRTimes, discount := currd]
    }
    return(nGramTable)
}
oneGramTable = goodTurningDiscountCoef(oneGramTable, k = 5)
twoGramTable = goodTurningDiscountCoef(twoGramTable, k = 5)
threeGramTable = goodTurningDiscountCoef(threeGramTable, k = 5)
fourGramTable = goodTurningDiscountCoef(fourGramTable, k = 3)

twoGramTable_leftOverProb = 
    twoGramTable[, .(leftoverprob =
                         calcLeftOverProb(freq, discount)), by=firstTerms]
threeGramTable_leftOverProb = 
    threeGramTable[, .(leftoverprob =
                           calcLeftOverProb(freq, discount)), by=firstTerms]
fourGramTable_leftOverProb = 
    fourGramTable[, .(leftoverprob =
                          calcLeftOverProb(freq, discount)), by=firstTerms]

write.csv(oneGramTable, 'oneGramTable_Ratz.csv', row.names = F)
write.csv(twoGramTable, 'twoGramTable_Ratz.csv', row.names = F)
write.csv(threeGramTable, 'threeGramTable_Ratz.csv', row.names = F)
write.csv(fourGramTable, 'fourGramTable_Ratz.csv', row.names = F)
########################################################################
#Start here if nothing needs to be changed
oneGramTable = fread("oneGramTable.csv")
twoGramTable = fread("twoGramTable.csv")
threeGramTable = fread("threeGramTable.csv")
fourGramTable = fread("fourGramTable.csv")

