library(tm)
library(SnowballC)
library(slam)
library(readr)
library(data.table)
library(dplyr)
#start here when model is fixed
oneGramTable = fread("oneGramTable.csv")
twoGramTable = fread("twoGramTable.csv")
threeGramTable = fread("threeGramTable.csv")
#fourGramTable = fread("fourGramTable.csv")
#if use smaller dictionary the algrithm will run much faster
#trim_oneGramTable = oneGramTable[freq > 1]

# function to get prefix and tail word
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

#  In a 3gram model, get probability of occurance of all the possible last term (observed and unobserved)
#  unobserved means not observed in 3gram dictionary, but may exist in 2grams or 1grams dictironary


#______________________________________________________________________
#Step 1 Calculate Probabilities of Words Completing Unobserved Trigrams
## Returns a two column data.frame of observed trigrams that start with the
## bigram prefix in the first column named firstTerms and
## frequencies/counts in the second column named freq. If no observed trigrams
## that start with firstTerms exist, an empty data.frame is returned.

#calculate the conditional prob taking account of the discount factor (gama3 = 0.5)
getObs3gramProb <- function(obs3grams, twoGramTable, prefix3gram, triDisc=0.5) {
    if(nrow(obs3grams) < 1) return(NULL)
    prefix3gramCount <-twoGramTable[firstTerms == getFirstTerms(prefix3gram,1) 
                                  & lastTerm == getLastTerms(prefix3gram,1)]$freq
    obs3grams$prob <- log10((obs3grams$freq - triDisc) / prefix3gramCount)
    return(obs3grams)
}


#############################################################################################
#Step 2 Calculate Probabilities of Words Completing Unobserved Trigrams
#___________________________________________________________________________________________
#i) Find all the words that complete unobserved trigrams. 
#  These are the words in the set for wi:(wi−2,wi−1)wescribed earlier.

#getUnobs3gramTails takes the generated observed 3gram data frame, and the oneGramTable(as dictionary),
#and output the words that do not exist in the 3gram tail words
getUnobs3gramTails <- function(Obs3grams, oneGramTable) {
    obs_3gram_tails <- Obs3grams$lastTerm
    unobs_3gram_tails <- oneGramTable[!(oneGramTable$firstTerms %in% obs_3gram_tails)]$firstTerms
    return(unobs_3gram_tails)
}


#____________________________________________________________________________________________
#ii). Calculate discounted (left-over) probability mass at the bigram level α(wi−1)
getAlpha2gram <- function(prefix2gram, oneGramTable, twoGramTable, bigDisc=0.5) {
    # get all bigrams that start with wi-1
    # prefix2gram (i.e. wi-1) is the second word in the prefix3gram (the word you get
    # after tropping off the oldest word (i.e. wi-2) from prefix3gram)
    prefix2gramCount = oneGramTable[firstTerms == prefix2gram]$freq 
    oneGroupIn2Gram = twoGramTable[firstTerms == prefix2gram]
    if(nrow(oneGroupIn2Gram) < 1) {
 #       print('wi-1 not found') 
        return(1)
          }
    # if no bigram starts with wi-1, return alpha = 1!!!
    alpha2gram <- 1 - sum(oneGroupIn2Gram$freq - bigDisc) / prefix2gramCount
    return(alpha2gram)
}

#_____________________________________________________________________________________________
#iii). back-off prob for observed and unobserved 2grams.
# column prefix2gram (w_i-1) and colum unobs3gramTails form a dataframe (bo2grams), it contains all the  
# bigrams for which we want to calculate back-off prob

getObs2gramBoProbs <- function(obsBo2grams, oneGramTable, bigDisc=0.5) {
    if(nrow(obsBo2grams)<1) return(NULL)
    prefix2gramCount <-oneGramTable[firstTerms %in% obsBo2grams$firstTerms]$freq
    obsBo2grams$prob <- log10((obsBo2grams$freq - bigDisc) / prefix2gramCount)
    return(obsBo2grams)
}

## Returns a dataframe of 2 columns: ngram and prob.  Values in the ngram
## column are unobserved bigrams of the form: w2_w1.  The values in the prob
## column are the backed off probability estimates q_bo(w1 | w2) calculated
## from from equation 16.

getUnobs2gramBoProbs<- function(unobBo2grams, oneGramTable, alpha2gram) {
    qboUnobs2gram <- oneGramTable[firstTerms %in% unobBo2grams$lastTerm]
    alfreq <- sum(qboUnobs2gram$freq)
    qboUnobs2gram$prob <- log10(alpha2gram * qboUnobs2gram$freq / alfreq)
    return(qboUnobs2gram)
}

#just check that the total prbo of unobserved bigrams should add up to leftover prob
#alpha_2gram (ie.g. sum(10^Q_UnobsBo2grams$prob) == alpha_2gram)
#___________________________________________________________________________________________
#iv) calculate leftover probability at 3gram level: alpha_3gram, α(wi−2,wi−1)

## Returns the total probability mass discounted from all observed trigrams.
## calculated from equation 14. This is the amount of probability mass which is
## redistributed to UNOBSERVED trigrams. If no trigrams starting with
## bigram$ngram[1] exist, 1 is returned.

getAlpha3gram <- function(prefix3gram, twoGramTable, threeGramTable, triDisc=0.5) {
    prefix3gramCount <-twoGramTable[firstTerms == getFirstTerms(prefix3gram,1) 
                                    & lastTerm == getLastTerms(prefix3gram,1)]$freq
    oneGroupIn3Gram = threeGramTable[firstTerms == prefix3gram]
    if(nrow(oneGroupIn3Gram) < 1) { 
#        print('wi-2:i-1 not found')
        return(1)
        }
    # if no trigram starts with prefix3gram, return alpha = 1
    alpha3gram <- 1 - (sum(oneGroupIn3Gram$freq - triDisc) / prefix3gramCount)
    return(alpha3gram)
}

#___________________________________________________________________________________________
#v) calculate the unobserved 3gram backoff probability
getUnobs3gramProbs <- function(prefix3gram, qbo2grams, alpha3gram) {
    sum_2gram_Probs = sum(10^qbo2grams$prob)
    qbo2grams$prob <- log10(alpha3gram * (10^qbo2grams$prob) / sum_2gram_Probs)
    return(qbo2grams)
}

#______________________________________________________________________________________________
# Step 3: finally get the predict function
getPrediction <- function(qbo3grams) {
    # pull off tail word of highest prob trigram
    first_maxProb = max(qbo3grams$prob)
    second_maxProb = max(qbo3grams[prob != first_maxProb]$prob)
    third_maxProb = max(qbo3grams[(prob != first_maxProb)&(prob != second_maxProb)]$prob)
    maxProb = c(first_maxProb, second_maxProb, third_maxProb)
    prediction = data.table(Prediction = qbo3grams[prob %in% maxProb]$lastTerm, 
                            Probability = 10^qbo3grams[prob %in% maxProb]$prob)
#    print('Are you think of...')
    return(prediction)
}
## implement step 1-3
## get observed 3gram backoff prob
## bigger gamma2 and gamm3 means leftover probability is bigger for lower order grams (more weights on
## lower order grams)
## gamma2 = c(0.2, 0.4, 0.6, 0.8); gamma3 = c(0.2, 0.4, 0.6, 0.8) could be optimzed by CV
predict3gramTailword = function(input, gamma2=0.5, gamma3=0.5) {
    prefix3gram = getNgrams(input, 2)
    obs_3grams = threeGramTable[firstTerms == prefix3gram]
    obs_3gram_prob = getObs3gramProb(obs_3grams, twoGramTable, prefix3gram, gamma3)
    ## get unobserved 3gram tail
    unobs_3gram_tails <- getUnobs3gramTails(obs_3grams, oneGramTable)
    ## within the unobserved 3grams, after tropping of the oldest word, get observed and unobserved 2grams
    ## and get bockoff prob for all the possible 2grams
    prefix2gram <- getLastTerms(prefix3gram, 1)  #trop off the oldest word (i.e. wi-2)
    alpha_2gram <- getAlpha2gram(prefix2gram, oneGramTable, twoGramTable, gamma2)
    bo2grams = data.table(firstTerms = prefix2gram, lastTerm = unobs_3gram_tails)
    obsBo2grams = twoGramTable[(firstTerms %in% bo2grams$firstTerms)
                               &(lastTerm %in% bo2grams$lastTerm)]
    unobsBo2grams = bo2grams[!(lastTerm%in%obsBo2grams$lastTerm)]
    Q_obsBo2grams = getObs2gramBoProbs(obsBo2grams, oneGramTable)
    Q_UnobsBo2grams = getUnobs2gramBoProbs(unobsBo2grams, oneGramTable, alpha_2gram)
    Q_UnobsBo2grams$lastTerm = Q_UnobsBo2grams$firstTerms
    Q_UnobsBo2grams$firstTerms = prefix2gram
    qbo_2grams <- rbind(Q_obsBo2grams, Q_UnobsBo2grams)
    qbo_2grams = qbo_2grams[order(qbo_2grams$prob, decreasing = TRUE),]
    ## calculate alpha at 3gram level
    alpha_3gram <- getAlpha3gram(prefix3gram, twoGramTable, threeGramTable, gamma3)
    unobs_3gram_prob = getUnobs3gramProbs(prefix3gram, qbo_2grams, alpha_3gram)
    ## generate dataframe storing probability for observed and unobserved 3grams
    qbo_3grams <- rbind(obs_3gram_prob, unobs_3gram_prob)
    qbo_3grams <- qbo_3grams[order(-qbo_3grams$prob), ]  
    ## prediction is the 3gram tail with highest prob
    getPrediction(qbo_3grams)
}

getProb3gramModel = function(triGram, gamma2=0.5, gamma3=0.5) {
    prefix3gram = getFirstTerms(triGram, n=2)
    prefix2gram = getLastTerms(prefix3gram, 1)  #trop off the oldest word (i.e. wi-2)
    tail3gram = getLastTerms(triGram, n=1)
    alpha_2gram <- getAlpha2gram(prefix2gram, oneGramTable, twoGramTable, gamma2) 
    #alpha_2gram=1 means wi-1 not found
    alpha_3gram <- getAlpha3gram(prefix3gram, twoGramTable, threeGramTable, gamma3)
    #alpha_3gram=1 means wi-2:i-1 not found
    
    obs_3grams = threeGramTable[firstTerms == prefix3gram]
    ## most simple senario: found observed 3grams
    if (nrow(obs_3grams)>0) {
        obs_3gram_prob = getObs3gramProb(obs_3grams, twoGramTable, prefix3gram, gamma3)
        record_obs_3gram =  obs_3gram_prob[firstTerms == prefix3gram & lastTerm == tail3gram]
        if (nrow(record_obs_3gram)>0) {
            #            print ('observed tri-gram found')
            return(10^record_obs_3gram$prob) } 
    } 
    unobs_3gram_tails <- getUnobs3gramTails(obs_3grams, oneGramTable)
    bo2grams = data.table(firstTerms = prefix2gram, lastTerm = unobs_3gram_tails)
    obsBo2grams = twoGramTable[(firstTerms %in% bo2grams$firstTerms)
                               &(lastTerm %in% bo2grams$lastTerm)]
    unobsBo2grams = bo2grams[!(lastTerm%in%obsBo2grams$lastTerm)]
    Q_obsBo2grams = getObs2gramBoProbs(obsBo2grams, oneGramTable)
    Q_UnobsBo2grams = getUnobs2gramBoProbs(unobsBo2grams, oneGramTable, alpha_2gram)
    Q_UnobsBo2grams$lastTerm = Q_UnobsBo2grams$firstTerms
    Q_UnobsBo2grams$firstTerms = prefix2gram
    qbo_2grams <- rbind(Q_obsBo2grams, Q_UnobsBo2grams)
    qbo_2grams = qbo_2grams[order(qbo_2grams$prob, decreasing = TRUE),]
    unobs_3gram_prob = getUnobs3gramProbs(prefix3gram, qbo_2grams, alpha_3gram)
    record_unobs_3gram = unobs_3gram_prob[(firstTerms == prefix2gram & lastTerm == tail3gram)]
    if (nrow(record_unobs_3gram) > 0) {
        #        print ('unobserved tri-gram probability estimated as') 
        return(10^record_unobs_3gram$prob) 
    } else {
        print ('no uni-gram found!')
        return (10^min(unobs_3gram_prob$prob))
    }
}

# get n-gram from an input sentence, if n<=length(sentence), get the n-gram, if n>length(sentence) get
# all the words and put a starting symbol at the begining
gramCount = function (inputString) {
    # inputString = gsub('[[:punct:]]', "", inputString)
    # inputString = gsub('[0-9]', "", inputString)
    wordList = unlist(strsplit(inputString, " ", fix = T))
    #wordList = wordList[wordList!='']
    length(wordList)
}

getNgrams = function (inputString, n=2) {
    # get the last ngram from a input sentence
    inputString = gsub('[[:punct:]]', "", inputString)
    inputString = gsub('[0-9]', "", inputString)
    wordList = unlist(strsplit(inputString, " ", fix = T))
    wordList = wordList[wordList!=''] #strip extra white space
    l = length(wordList)
    if (l>=n) ngram = paste(c(wordList[(l-n+1):l]), collapse = ' ')
    else {
        ngram = paste(c(wordList[1:l]), collapse = ' ')
        ngram = paste('<s>', ngram, collapse = ' ')
    }
    tolower(ngram)
}

