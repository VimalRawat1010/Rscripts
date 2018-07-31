This is the code we use in this section. Just copy and paste it into R and follow along with the videos.

#####################################

# Social Media Mining

# Twitter is a great source for sentiment data and social media mining
# furthermore it is quite easy to get significant amounts of data
# to be able to scrape data from Twitter you need a standard Twitter account
# and you need to update it to a developer account
# note that Twitter limits the amount of searches you can perform (15min: 15 scrapes)

# package twitteR
library("twitteR")


install.packages(c('ROAuth','RCurl'))
require('ROAuth')
require('RCurl')

library(RCurl)
consumer_key <- "tD2bCFAXCtdHoFABxCXW5MSEI"
consumer_secret <- "H5LQ7UbjTvRp9HJKEU3ixzaKbV38OrI7zdCfPqhbSra2a6Uswl"
access_token <- "1103421396-yyE408tG2g9Ce5WaO3Mr7q9BAsBfX9Za6sEBMBX"
access_secret <- "OYFmPUMrmJMbLfDZSjPdkulE2VOmZ6Flcn8dk871hrxvj"

setup_twitter_oauth(consumer_key, consumer_secret, access_token, access_secret)




library("twitteR")

# we need to specify the cainfo to avoid a SSL cert error - this is for Windows machines
# Lets check the latest tweets of Udemy
userTimeline("Udemy", cainfo="cacert.pem")

# searchTwitter is the main function of the package

?searchTwitter

# arguments: since and until are for time specifications
# lang: for languge specification
# geocode: for location specification

# we are now scraping 1k tweekts for Udemy, and we als specify our certificate
udemytweets = searchTwitter("#Udemy", n=1000)

# as you can see, scraping that data is quite time consuming - your machine limits the
# the efficiency and speed of your mining 
#
