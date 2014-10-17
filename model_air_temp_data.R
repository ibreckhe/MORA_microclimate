## Script to fill gaps in air temperature data.
## Author: Ian Breckheimer
## Date: 16 October 2014

#### Sets up workspace and load data ####

## Loads required packages.
library(spBayes)
library(MBA)
library(geoR)

## Loads and formats data.
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
topowx <- read.csv("topowx_daily_2004_2012.csv")
topowx$DATE <- as.Date(topowx$DATE)

prism <- read.csv("PRISM_daily_2004_2014.csv")
prism$DATE <- as.Date(prism$DATE)


setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/cleaned_dailyair/")
tavg <- read.csv("alldat_daily_tavg.csv")
tavg$DATE <- as.Date(tavg$DATE)
tmax <- read.csv("alldat_daily_tmax.csv")
tmax$DATE <- as.Date(tmax$DATE)
tmin <- read.csv("alldat_daily_tmin.csv")
tmin$DATE <- as.Date(tmin$DATE)
meta <- read.csv("site_metadata.csv")
meta$alt_code <- gsub(pattern="-",replacement=".",x=meta$location,fixed=TRUE)

#### Data exploration: what does the pattern of spatial dependence look like? ####

##Gets tavg data for September 24th, 2013.
date_range <- 2600:2800

## Vectors to hold linear model coefficients and residual spatial parameters.
intercepts <- rep(NA,length(date_range))
lapses <- rep(NA,length(date_range))
canopies <- rep(NA,length(date_range))
relevs <- rep(NA,length(date_range))
rsqs <- rep(NA,length(date_range))
aics <- rep(NA,length(date_range))
maes <- rep(NA,length(date_range))
nuggets <- rep(NA,length(date_range))
ranges <- rep(NA,length(date_range))

for (i in date_range){
  index <- (i - min(date_range))+1
  day <- t(tmax[i,2:dim(tmax)[2]])
  colnames(day) <- "temp"
  day_meta <- merge(data.frame(alt_code=rownames(day),tmax=day),meta)
  #plot(temp~elev,data=day_meta,main=tavg$DATE[date_range[index]])
  #text(x=day_meta$elev,y=day_meta$temp,labels=day_meta$alt_code)
  test_lm <- lm(temp~elev+I(asin(ccov_81m))+relev_729m+utmx+utmy,data=day_meta)
  intercepts[index] <- test_lm$coefficients[1]
  lapses[index] <- test_lm$coefficients[2]
  canopies[index] <- test_lm$coefficients[3]
  relevs[index] <- test_lm$coefficients[4]
  rsqs[index] <- summary(test_lm)$r.squared
  aics[index] <- AIC(test_lm)
  summary(test_lm)
  temp_resid <- data.frame(day_meta[complete.cases(day_meta),],residuals(test_lm))
  maes[index] <- mean(abs(temp_resid$residuals.test_lm.))
  coordinates(temp_resid) <- ~utmx+utmy
  
  ##Residual variogram.
  resid_variog <- variog(coords=coordinates(temp_resid),data=temp_resid$residuals.test_lm.,
                         max.dist=20000,option="bin",breaks=seq(0,20000,by=1000),messages=FALSE)
  resid_fit <- variofit(resid_variog,cov.model="spherical",fix.nugget=FALSE,max.dist=20000,messages=FALSE)
  nuggets[index] <- resid_fit$cov.pars[1]
  ranges[index] <- resid_fit$cov.pars[2]
  #plot(resid_variog)
  #lines(resid_fit)
}

pdf("../results/lm_coefficients_diagnostics_tmax.pdf",width=10,height=10)
par(mfrow=c(8,1),mar=c(0,4,1,1),oma=c(4,0,0,0))
plot(tavg$DATE[date_range],type="l",xaxt="n",lapses,ylab="Lapse Rate",ylim=c(-0.01,0.01))
abline(h=0,lty=2)
text(x=as.numeric(max(tavg$DATE[date_range])),y=7,labels=round(mean(lapses),digits=3))
plot(tavg$DATE[date_range],type="l",xaxt="n",canopies,ylab="Can. Cov.",ylim=c(-8,8))
abline(h=0,lty=2)
text(x=as.numeric(max(tavg$DATE[date_range])),y=3,labels=round(mean(canopies),digits=3))
plot(tavg$DATE[date_range],type="l",xaxt="n",relevs,ylab="R. elev.",ylim=c(-0.1,0.1))
abline(h=0,lty=2)
text(x=as.numeric(max(tavg$DATE[date_range])),y=0.08,labels=round(mean(relevs),digits=3))
plot(tavg$DATE[date_range],type="l",xaxt="n",rsqs,ylim=c(0,1),ylab=expression(r^2))
text(x=as.numeric(max(tavg$DATE[date_range])),y=0.1,labels=round(mean(rsqs),digits=3))
plot(tavg$DATE[date_range],type="l",xaxt="n",aics,ylim=c(0,400),ylab="AIC")
text(x=as.numeric(max(tavg$DATE[date_range])),y=50,labels=round(mean(aics),digits=3))
plot(tavg$DATE[date_range],type="l",xaxt="n",maes,ylab="Mean absolute error (C)",ylim=c(0,2))
text(x=as.numeric(max(tavg$DATE[date_range])),y=0.25,labels=round(mean(maes),digits=3))
plot(tavg$DATE[date_range],type="l",xaxt="n",nuggets,ylab="Residual Nugget")
text(x=as.numeric(max(tavg$DATE[date_range])),y=2,labels=round(mean(nuggets),digits=3))
plot(tavg$DATE[date_range],type="l",ranges,ylab="Residual Range")
dev.off()

#### Initial Models: Spatial model fitting: what are good starting values for the parameters? ####

#### Temporal Subsets: Fitting space-time model on a temporal subset of the data.####

#### Computation: fitting the space-time model on the full dataset.####

#### Prediction: Using the posterior estimates to predict data at new locations.####
