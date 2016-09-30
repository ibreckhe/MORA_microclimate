##Script to predict minimum soil temp maps.

library(lme4)
library(spBayes)
library(ggplot2)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/processed/")

##Reads in data.
snowdat <- read.csv("microclimate_snow_coords_final_10_19_15.csv")

ggplot(snowdat)+
  geom_point(aes(x=elevation,y=min_soil_temp,fill=water_year))+
  theme_bw()