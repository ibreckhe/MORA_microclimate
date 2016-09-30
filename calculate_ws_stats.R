##Script to calculate disparity, sensitivity, and decoupling for weather station data.
##Author: Ian Breckheimer
##Date: May 3rd, 2016

library(rnoaa)
library(plyr)
library(reshape2)
options(noaakey = "NjuEQnGAERfxLrZNIiHqMUJsSeqiPuRR")

##Load weather station data.
ws <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/GHCN/GHCN_rainier_all_1981_2015.csv")
ws$date <- as.Date(as.character(ws$DATE),format="%Y%m%d")
ws_meta <- unique(data.frame(ws$STATION,ws$STATION_NAME,ws$LATITUDE,ws$LONGITUDE,ws$ELEVATION))
write.csv(ws_meta,"~/Code/mora_climate/GHCN_locations.csv")
