# This function is used to get the probability of a given text, using Katz Backoff 
# (with Good-Turing Discounting).

getNgrams = function (inputSentence, n) {
    # get the last ngram from a input sentence
    inputSentence = gsub('[[:punct:]]', "", inputSentence)
    wordList = unlist(strsplit(inputSentence, " ", fix = T))
    l = length(wordList)
    if (l>=n) ngram = paste(c(wordList[(l-n+1):l]), collapse = ' ')
    
    else {ngram = paste(c(wordList[1:l]), collapse = ' ')
          ngram = paste('<s>', ngram, collapse = ' ')
              }
    return(tolower(ngram))
}
###################################################
getProbabilityFrom2Gram = function(inputString){
    # Preprocessing
    input = getNgrams(inputString, 2)
    inFirstTerms2gram = getFirstTerms(input, 1)
    inLastTerm2gram = getLastTerms(input, 1)
    #initializing
    finalProb = 0
    #locate 3-grams where firstTerms match
    oneGroupIn2Gram = twoGramTable[firstTerms == inFirstTerms2gram]
    if (nrow(oneGroupIn2Gram) > 0){
        # Algorithm here
        oneRecordIn2Gram = twoGramTable[firstTerms == inFirstTerms2gram 
                                          & lastTerm == inLastTerm2gram]
        if (nrow(oneRecordIn2Gram) > 0){
            # We found one in 3-gram
            all_freq = sum(oneGroupIn2Gram$freq)
            finalProb = ((oneRecordIn2Gram$discount * oneRecordIn2Gram$freq) / all_freq)
            finalProb = log10(finalProb+1)
            print('term found in 2-gram table')
            ### We're done!
        } else {
            # We only have hope in 1-gram!
            # # Get the left-over prob. so that we can distribute it for lower-order grams.
            beta_leftoverprob = 
                 twoGramTable_leftOverProb[firstTerms == inFirstTerms2gram]$leftoverprob
            oneGroupIn1Gram = oneGramTable  # we don't have "firstTerms" here!
            oneRecordIn1Gram = oneGramTable[firstTerms == inLastTerm2gram] 
            oneGroupIn1Gram_Remain =
                oneGroupIn1Gram[(oneGroupIn1Gram$firstTerms%in%oneGroupIn2Gram$lastTerm)]
            
            all_freq = sum(oneGroupIn1Gram$freq)
            alpha = beta_leftoverprob/(1-
                sum((oneGroupIn1Gram_Remain$freq*oneGroupIn1Gram_Remain$discount)/
                        all_freq))
            finalProb = alpha*((oneRecordIn1Gram$freq*oneRecordIn1Gram$discount)/ all_freq)
            finalProb = log10(finalProb+1)
            print('term found in uni-gram table')
            
        }
     
    } else {
        sprintf("[%s] not found in the 2-gram model.", inFirstTerms2gram)
        finalProb = 0
    }
    return(finalProb)
}
#########################################################################################

getProbabilityFrom3Gram = function(inputString){
    # Preprocessing
    input = getNgrams(inputString, 3)
    inFirstTerms3gram = getFirstTerms(input, 2)
    inLastTerm3gram = getLastTerms(input, 1)
    #initializing
    finalProb = 0
    #locate 3-grams where firstTerms match
    oneGroupIn3Gram = threeGramTable[firstTerms == inFirstTerms3gram]
    if (nrow(oneGroupIn3Gram) > 0){
        # Algorithm here
        oneRecordIn3Gram = threeGramTable[firstTerms == inFirstTerms3gram 
                                          & lastTerm == inLastTerm3gram]
        if (nrow(oneRecordIn3Gram) > 0){
            # We found one in 3-gram
            all_freq = sum(oneGroupIn3Gram$freq)
            finalProb = ((oneRecordIn3Gram$discount * oneRecordIn3Gram$freq) / all_freq)
            finalProb = log10(finalProb+1)
            print('term found in 3-gram table')
            ### We're done!
        } else {
            # NOT found in 3-gram => check 2-gram & 1-gram
            input = getNgrams(inputString, n = 2)
            inFirstTerms2gram = getFirstTerms(input, 1)
            inLastTerm2gram = getLastTerms(input, 1)
            
            # Get the left-over prob. so that we can distribute it for lower-order grams.
            beta_leftoverprob = 
                threeGramTable_leftOverProb[firstTerms == inFirstTerms3gram]$leftoverprob
            
            oneGroupIn2Gram = twoGramTable[firstTerms == inFirstTerms2gram]
            oneRecordIn2Gram = twoGramTable[firstTerms == inFirstTerms2gram 
                                            & lastTerm == inLastTerm2gram]
            if (nrow(oneRecordIn2Gram) > 0){
                # We found one in 2-gram!
                # We only consider ones that do not appear in 3-grams...
                oneGroupIn2Gram_Remain = 
                    oneGroupIn2Gram[(oneGroupIn2Gram$lastTerm%in%oneGroupIn3Gram$lastTerm)]
                
                all_freq = sum(oneGroupIn2Gram$freq)
                #alpha is back-off weight
                alpha = beta_leftoverprob/(1-
                        sum((oneGroupIn2Gram_Remain$freq*oneGroupIn2Gram_Remain$discount)
                        / all_freq))
                
                finalProb = alpha*((oneRecordIn2Gram$freq * oneRecordIn2Gram$discount) /
                            all_freq) #done!
                finalProb = log10(finalProb+1)
                print('term found in 2-gram table')
            } else {
                # # Get the left-over prob. so that we can distribute it for lower-order grams.
                # beta_leftoverprob = 
                #     twoGramTable_leftOverProb[firstTerms == inFirstTerms2gram]$leftoverprob
                
                oneGroupIn1Gram = oneGramTable # we don't have "firstTerms" here!
                oneRecordIn1Gram = oneGramTable[firstTerms == inLastTerm2gram] 
                # what if this returns "zero" row?
                oneGroupIn1Gram_Remain =
                    oneGroupIn1Gram[(oneGroupIn1Gram$firstTerms%in%oneGroupIn3Gram$lastTerm)]
                all_freq = sum(oneGroupIn1Gram$freq)
                
                alpha = beta_leftoverprob/(1-
                        sum((oneGroupIn1Gram_Remain$freq*oneGroupIn1Gram_Remain$discount)/
                            all_freq))
                
                finalProb = alpha*((oneRecordIn1Gram$freq*oneRecordIn1Gram$discount)
                            / all_freq)
                finalProb = log10(finalProb+1)
                print('term found in uni-gram table')
            }
        }
    } else {
        sprintf("[%s] not found in the 3-gram model.", inFirstTerms3gram)
        finalProb = 0
    }
    return(finalProb)
    
}
#########################################################################################
getProbabilityFrom4Gram = function(inputString){
    input = getNgrams(inputString, 4)
    inFirstTerms4gram = getFirstTerms(input, 3)
    inLastTerm4gram = getLastTerms(input, 1)
    #initializing
    finalProb = 0
    #locate 4-grams where firstTerms match
    oneGroupIn4Gram = fourGramTable[firstTerms == inFirstTerms4gram]
    if (nrow(oneGroupIn4Gram) > 0){
        oneRecordIn4Gram = fourGramTable[firstTerms == inFirstTerms4gram 
                                          & lastTerm == inLastTerm4gram]
        if (nrow(oneRecordIn4Gram) > 0){
            # We found one in 4-gram
            all_freq = sum(oneGroupIn4Gram$freq)
            finalProb = ((oneRecordIn4Gram$discount * oneRecordIn4Gram$freq) / all_freq)
            finalProb = log10(finalProb+1)
            print("term found in the 4-gram table.")
        } else {
            # NOT found in 4-gram => check 3, 2, 1-gram 
            input = getNgrams(inputString, n = 3)
            inFirstTerms3gram = getFirstTerms(input, 2)
            inLastTerm3gram = getLastTerms(input, 1)
            beta_leftoverprob = 
                fourGramTable_leftOverProb[firstTerms == inFirstTerms4gram]$leftoverprob
            oneGroupIn3Gram = threeGramTable[firstTerms == inFirstTerms3gram]
            oneRecordIn3Gram = threeGramTable[firstTerms == inFirstTerms3gram 
                                            & lastTerm == inLastTerm3gram]
            if (nrow(oneRecordIn3Gram) > 0){
                # We found one in 3-gram!
                # We only consider ones that do not appear in 4-grams...
                oneGroupIn3Gram_Remain = 
                    oneGroupIn3Gram[(oneGroupIn3Gram$lastTerm%in%oneGroupIn4Gram$lastTerm)]
                all_freq = sum(oneGroupIn3Gram$freq)
                alpha = beta_leftoverprob/(1-
                        sum((oneGroupIn3Gram_Remain$freq*oneGroupIn3Gram_Remain$discount)
                        / all_freq))
                
                finalProb = alpha*((oneRecordIn3Gram$freq * oneRecordIn3Gram$discount) /
                                       all_freq) #done!
                finalProb = log10(finalProb+1)
                print("term found in the 3-gram table.")
            } else {
                # NOT found in 3-gram => check 2, 1-gram
                input = getNgrams(inputString, n = 2)
                inFirstTerms2gram = getFirstTerms(input, 1)
                inLastTerm2gram = getLastTerms(input, 1)
                oneGroupIn2Gram = twoGramTable[firstTerms == inFirstTerms2gram]
                oneRecordIn2Gram = twoGramTable[firstTerms == inFirstTerms2gram 
                                                & lastTerm == inLastTerm2gram]
                if (nrow(oneRecordIn2Gram) > 0){
                    oneGroupIn2Gram_Remain = 
                        oneGroupIn2Gram[(oneGroupIn2Gram$lastTerm%in%oneGroupIn4Gram$lastTerm)]
                    all_freq = sum(oneGroupIn2Gram$freq)
                    alpha = beta_leftoverprob/(1-
                        sum((oneGroupIn2Gram_Remain$freq*oneGroupIn2Gram_Remain$discount)
                            / all_freq))
                    
                    finalProb = alpha*((oneRecordIn2Gram$freq * oneRecordIn2Gram$discount) /
                                           all_freq)
                    finalProb = log10(finalProb+1)
                    print("term found in the 2-gram table.")
                } else {
                    oneGroupIn1Gram = oneGramTable 
                    oneRecordIn1Gram = oneGramTable[firstTerms == inLastTerm2gram] 
                    oneGroupIn1Gram_Remain =
                        oneGroupIn1Gram[(oneGroupIn1Gram$firstTerm%in%oneGroupIn4Gram$lastTerm)]
                    all_freq = sum(oneGroupIn1Gram$freq)
                    
                    alpha = beta_leftoverprob/(1-
                        sum((oneGroupIn1Gram_Remain$freq*oneGroupIn1Gram_Remain$discount)/
                                all_freq))
                    
                    finalProb = alpha*((oneRecordIn1Gram$freq*oneRecordIn1Gram$discount)
                                       / all_freq)
                    finalProb = log10(finalProb+1)
                    print("term found in the uni-gram table.")
                }
            }      
        }
        
    } else {
        sprintf("[%s] not found in the 4-gram model.", inFirstTerms4gram)
        finalProb = 0
    }
    return(finalProb)
    
}
    