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
oneGramTable=oneGramTable[with(oneGramTable, order(firstTerms, -freq)), ]

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
