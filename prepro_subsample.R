library(caret)
library(tm)
library(SnowballC)
library(slam)
library(ggplot2)
library(readr)

con<-file("en_US/en_US.twitter.txt", encoding = 'UTF-8')
twitter<-read_lines(con)
con<-file("en_US/en_US.news.txt", encoding = 'UTF-8')
news<-read_lines(con)
con<-file("en_US/en_US.blogs.txt", encoding = 'UTF-8')
blogs<-read_lines(con)

#a function that remove examples with dirty words and non-english words 
sampleCleaning<-function(document) {
    #remove dirty words
    document<-gsub("[Ff][Uu][Cc][Kk]", "", document, ignore.case = TRUE)
    
    #remove non-english words
    document<-gsub("dat2", "", iconv(document, "latin1", "ASCII", sub="dat2"))
    
}
#clean and sub_sample
twitter1 = sampleCleaning(twitter)
news1 = sampleCleaning(news)
blogs1 = sampleCleaning(blogs)
write(twitter1, "twitter.txt")
write(news1, "news.txt")
write(blogs1, "blogs.txt")

subIdx_twitter<-rbinom(n = length(twitter1), size = 1, prob = 0.05) 
subIdx_news<-rbinom(n = length(news1), size = 1, prob = 0.1)
subIdx_blogs<-rbinom(n = length(blogs), size = 1, prob = 0.1)
sub_twitter<-twitter1[which(subIdx_twitter!=0)]
sub_news<-news1[which(subIdx_news!=0)]
sub_blogs<-blogs1[which(subIdx_blogs!=0)]

sample = c(sub_twitter, sub_news, sub_blogs)
write(sample, "sub_sample.txt")
write(sub_twitter, "sub_twitter.txt")
write(sub_news, "sub_news.txt")
write(sub_blogs, "sub_blogs.txt")

