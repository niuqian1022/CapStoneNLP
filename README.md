# CapStoneNLP
When you type a few words in a search engine or editing software, for example, being able to predict what's comming up next is a nice feature to have! The goal of this project is to predict next word occuring in text document based on previously seen words. 

Text samples were collected from Twitter, US news and on-line blogs. I implemented 3-gram models, and accouted for unseen n-grams using Ratz's back-off algrithm with absolute discounting. I built 'Word Scryer' which is a mini shinny app, and you are welcome to try out a few examples.
https://niuqian.shinyapps.io/word_scryer/

How to use the App.
Take a look at the tabs which include instruction, result and info for your reference.
Type a short sentence/phrase you want to begin with (not case sensitive).
Click Result and have fun.
A few second delay should be expected.

Caveats:
Discouting parameters should be optimized with cross validation;
Consider to remove some singleton from dictionary to speed up the algrithm;
Choices of higher order N-gram models versus simpler discouting methods versus larger sample size


