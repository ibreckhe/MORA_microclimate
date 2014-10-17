##Script to download Weather Underground weather station data for the vicinity of Mt. Rainier National Park
##Author: Ian Breckheimer
##Date: 30 September, 2014.

##sets up worspace
library(plyr)
library(reshape2)
library(jsonlite)
library(RCurl)

##Sets working directory.
setwd("~/code/MORA_microclimate")

##Constructs a URL to download station metadata for the vicinity of Mt. Rainier.
"http://api.wunderground.com/api/15198a882fc52abb/history_20140924/q/46.853207,-121.765361.json"