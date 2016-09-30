##Script to generate daily temperature anomaly time-series for all stations.
library(xts)
library(dplyr)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
source("~/code/MORA_microclimate/air_temp_functions.R")

## Loads and formats data.
tavg <- read.csv("../cleaned/tavg_daily_all_stations_2009_2015.csv")
tavg$DATE <- as.Date(tavg$DATE)
tmin <- read.csv("../cleaned/tmin_daily_all_stations_2009_2015.csv")
tmin$DATE <- as.Date(tavg$DATE)
tmax <- read.csv("../cleaned/tmax_daily_all_stations_2009_2015.csv")
tmax$DATE <- as.Date(tmax$DATE)

## Subsets to target date range.
tavg_start <- as.Date("2009-10-01")
tavg_end <- as.Date("2015-9-30")
tavg <- tavg[tavg$DATE %in% seq(tavg_start,tavg_end,by="day"),]

tmin_start <- as.Date("2009-10-01")
tmin_end <- as.Date("2015-9-30")
tmin <- tmin[tmin$DATE %in% seq(tmin_start,tmin_end,by="day"),]

tmax_start <- as.Date("2009-10-01")
tmax_end <- as.Date("2015-9-30")
tmax <- tmax[tmax$DATE %in% seq(tmax_start,tmax_end,by="day"),]

##Reads seasonal estimates of macro-micro relationships from disk.
indices_dry_tmax <- read.csv("../cleaned/meta_indices_dry_tmax_2009_2015.csv")
indices_wet_tmin <- read.csv("../cleaned/meta_indices_wet_tmin_2009_2015.csv")
indices_tavg <- read.csv("../cleaned/meta_indices_tavg_2009_2015.csv")

indices_DJF_tavg <- read.csv("../cleaned/indices_DJF_tavg_2009_2015.csv")
indices_DJF_tmin <- read.csv("../cleaned/indices_DJF_tmin_2009_2015.csv")
indices_DJF_tmax <- read.csv("../cleaned/indices_DJF_tmax_2009_2015.csv")

indices_MAM_tavg <- read.csv("../cleaned/indices_MAM_tavg_2009_2015.csv")
indices_MAM_tmin <- read.csv("../cleaned/indices_MAM_tmin_2009_2015.csv")
indices_MAM_tmax <- read.csv("../cleaned/indices_MAM_tmax_2009_2015.csv")

indices_JJA_tavg <- read.csv("../cleaned/indices_JJA_tavg_2009_2015.csv")
indices_JJA_tmin <- read.csv("../cleaned/indices_JJA_tmin_2009_2015.csv")
indices_JJA_tmax <- read.csv("../cleaned/indices_JJA_tmax_2009_2015.csv")

indices_SON_tavg <- read.csv("../cleaned/indices_SON_tavg_2009_2015.csv")
indices_SON_tmin <- read.csv("../cleaned/indices_SON_tmin_2009_2015.csv")
indices_SON_tmax <- read.csv("../cleaned/indices_SON_tmax_2009_2015.csv")

##Reads in and subsets PRISM data to same date range.
prism <- read.csv("PRISM_daily_1_2004_9_2015.csv")
prism$DATE <- as.Date(prism$DATE)

##Computes seasonal averages from the PRISM data.
prism$month <- as.factor(format(prism$DATE,format="%m"))
prism$seas <- NA
prism$seas[prism$month %in% c("12","01","02")] <- "DJF"
prism$seas[prism$month %in% c("03","04","05")] <- "MAM"
prism$seas[prism$month %in% c("06","07","08")] <- "JJA"
prism$seas[prism$month %in% c("09","10","11")] <- "SON"
prism$seas <- factor(prism$seas,levels=c("DJF","MAM","JJA","SON"))

##Creates merged data frames of seasonal estimates.
doy <- c(15,106,198,289)
disp_tavg <- t(data.frame(DJF=indices_DJF_tavg$disp,
                        MAM=indices_MAM_tavg$disp,
                        JJA=indices_JJA_tavg$disp,
                        SON=indices_SON_tavg$disp,
                        row.names=indices_MAM_tmin$site))
disp_tavg <- data.frame(cbind(disp_tavg,doy))

disp_tmax <- t(data.frame(DJF=indices_DJF_tmax$disp,
                        MAM=indices_MAM_tmax$disp,
                        JJA=indices_JJA_tmax$disp,
                        SON=indices_SON_tmax$disp,
                        row.names=indices_MAM_tmin$site))
disp_tmax <- data.frame(cbind(disp_tmax,doy))

disp_tmin <-t(data.frame(DJF=indices_DJF_tmin$disp,
                        MAM=indices_MAM_tmin$disp,
                        JJA=indices_JJA_tmin$disp,
                        SON=indices_SON_tmin$disp,
                        row.names=indices_MAM_tmin$site))
disp_tmin <- data.frame(cbind(disp_tmin,doy))

days <- data.frame(date=seq(as.Date("2008-07-16"),as.Date("2016-10-16"),by="day"))
days$doy <- format(days$date,format="%j")

days_disp_tavg <- merge(days,disp_tavg,by="doy",all.x=TRUE,sort=FALSE)
days_disp_tavg <- days_disp_tavg[order(days_disp_tavg$date),]

days_disp_tmax <- merge(days,disp_tmax,by="doy",all.x=TRUE,sort=FALSE)
days_disp_tmax <- days_disp_tmax[order(days_disp_tmax$date),]

days_disp_tmin <- merge(days,disp_tmin,by="doy",all.x=TRUE,sort=FALSE)
days_disp_tmin <- days_disp_tmin[order(days_disp_tmin$date),]

##Does the same for sensitivity.
##Creates merged data frames of seasonal estimates.
doy <- c(15,106,198,289)
sens_tavg <- t(data.frame(DJF=indices_DJF_tavg$sens,
                          MAM=indices_MAM_tavg$sens,
                          JJA=indices_JJA_tavg$sens,
                          SON=indices_SON_tavg$sens,
                          row.names=indices_MAM_tmin$site))
sens_tavg <- data.frame(cbind(sens_tavg,doy))

sens_tmax <- t(data.frame(DJF=indices_DJF_tmax$sens,
                          MAM=indices_MAM_tmax$sens,
                          JJA=indices_JJA_tmax$sens,
                          SON=indices_SON_tmax$sens,
                          row.names=indices_MAM_tmin$site))
sens_tmax <- data.frame(cbind(sens_tmax,doy))

sens_tmin <-t(data.frame(DJF=indices_DJF_tmin$sens,
                         MAM=indices_MAM_tmin$sens,
                         JJA=indices_JJA_tmin$sens,
                         SON=indices_SON_tmin$sens,
                         row.names=indices_MAM_tmin$site))
sens_tmin <- data.frame(cbind(sens_tmin,doy))

days <- data.frame(date=seq(as.Date("2008-07-16"),as.Date("2016-10-16"),by="day"))
days$doy <- format(days$date,format="%j")

days_sens_tavg <- merge(days,sens_tavg,by="doy",all.x=TRUE,sort=FALSE)
days_sens_tavg <- days_sens_tavg[order(days_sens_tavg$date),]

days_sens_tmax <- merge(days,sens_tmax,by="doy",all.x=TRUE,sort=FALSE)
days_sens_tmax <- days_sens_tmax[order(days_sens_tmax$date),]

days_sens_tmin <- merge(days,sens_tmin,by="doy",all.x=TRUE,sort=FALSE)
days_sens_tmin <- days_sens_tmin[order(days_sens_tmin$date),]

##Loops through each record and linearly interpolates the estimates.
days_disp_tavg_interp <- interpolateMicro(days_disp_tavg)
days_disp_tmin_interp <- interpolateMicro(days_disp_tmin)
days_disp_tmax_interp <- interpolateMicro(days_disp_tmax)

days_sens_tavg_interp <- interpolateMicro(days_sens_tavg)
days_sens_tmin_interp <- interpolateMicro(days_sens_tmin)
days_sens_tmax_interp <- interpolateMicro(days_sens_tmax)

##Trims the estimates to the right date range.
date_start <- as.Date("2009-10-01")
date_end <- as.Date("2015-09-30")
date_seq <- seq(date_start,date_end,by="day")

days_disp_tavg_interp <- days_disp_tavg_interp[days_disp_tavg_interp$date %in% date_seq,]
days_disp_tmin_interp <- days_disp_tmin_interp[days_disp_tmin_interp$date %in% date_seq,]
days_disp_tmax_interp <- days_disp_tmax_interp[days_disp_tmax_interp$date %in% date_seq,]

days_sens_tavg_interp <- days_sens_tavg_interp[days_sens_tavg_interp$date %in% date_seq,]
days_sens_tmin_interp <- days_sens_tmin_interp[days_sens_tmin_interp$date %in% date_seq,]
days_sens_tmax_interp <- days_sens_tmax_interp[days_sens_tmax_interp$date %in% date_seq,]

##Trims the PRISM time-series.
days_prism <- prism[prism$DATE %in% date_seq,]

##Creates interpolated seasonal expectations for the PRISM data.
prism_gr <- group_by(days_prism,seas)
prism_seas <- summarise(prism_gr,tavg_seas=mean(TEMP,na.rm=TRUE),
                                 tmax_seas=mean(TMAX,na.rm=TRUE),
                                 tmin_seas=mean(TMIN,na.rm=TRUE))
prism_seas$doy <- doy

days_prism_seas <- merge(days,prism_seas[,-1],by="doy",all.x=TRUE,sort=FALSE)
days_prism_seas <- days_prism_seas[order(days_prism_seas$date),]
days_prism_exp <- interpolateMicro(days_prism_seas)

##Trims the expectations.
days_prism_exp <- days_prism_exp[days_prism_exp$date %in% date_seq,]

##Populates a matrix of PRISM expectations.
prism_exp_tavg <- days_disp_tavg_interp
prism_exp_tavg[,3:dim(prism_exp_tavg)[2]] <- NA
for(i in 3:dim(prism_exp_tavg)[2]){
  prism_exp_tavg[,i] <- days_prism_exp$tavg_seas
}

prism_exp_tmin <- days_disp_tmin_interp
prism_exp_tmin[,3:dim(prism_exp_tmin)[2]] <- NA
for(i in 3:dim(prism_exp_tmin)[2]){
  prism_exp_tmin[,i] <- days_prism_exp$tmin_seas
}

prism_exp_tmax <- days_disp_tmax_interp
prism_exp_tmax[,3:dim(prism_exp_tmax)[2]] <- NA
for(i in 3:dim(prism_exp_tmax)[2]){
  prism_exp_tmax[,i] <- days_prism_exp$tmax_seas
}


##Populates a matrix of PRISM measurements.
prism_tavg <- days_disp_tavg_interp
prism_tavg[,3:dim(prism_tavg)[2]] <- NA
for(i in 3:dim(prism_tavg)[2]){
  prism_tavg[,i] <- days_prism$TEMP
}

prism_tmin <- days_disp_tmin_interp
prism_tmin[,3:dim(prism_tmin)[2]] <- NA
for(i in 3:dim(prism_tmin)[2]){
  prism_tmin[,i] <- days_prism$TMIN
}

prism_tmax <- days_disp_tmax_interp
prism_tmax[,3:dim(prism_tmax)[2]] <- NA
for(i in 3:dim(prism_tmax)[2]){
  prism_tmax[,i] <- days_prism$TMAX
}

##Creates matrices of centered prism data.
tavg_exp <- as.matrix(prism_tavg[,-c(1,2)]) - as.matrix(prism_exp_tavg[,-c(1,2)])
tmin_exp <- as.matrix(prism_tmin[,-c(1,2)]) - as.matrix(prism_exp_tmin[,-c(1,2)])
tmax_exp <- as.matrix(prism_tmax[,-c(1,2)]) - as.matrix(prism_exp_tmax[,-c(1,2)])

##Creates matrices of deviations from PRISM expected values.
tavg_dev <- tavg_exp * as.matrix(days_sens_tavg_interp[,-c(1,2)]) + as.matrix(days_disp_tavg_interp[,-c(1,2)])  
tmin_dev <- tmin_exp * as.matrix(days_sens_tmin_interp[,-c(1,2)]) + as.matrix(days_disp_tmin_interp[,-c(1,2)])  
tmax_dev <- tmax_exp * as.matrix(days_sens_tmax_interp[,-c(1,2)]) + as.matrix(days_disp_tmax_interp[,-c(1,2)])  

##Adds back the PRISM expected values.
tavg_exp_prism <- tavg_dev + as.matrix(prism_exp_tavg[,-c(1,2)])
tmin_exp_prism <- tmin_dev + as.matrix(prism_exp_tmin[,-c(1,2)])
tmax_exp_prism <- tmax_dev + as.matrix(prism_exp_tmax[,-c(1,2)])

##Subtracts expected values to produce anomalies.
tavg_anom <- as.matrix(tavg[,-1]) - tavg_exp_prism
tmax_anom <- as.matrix(tmax[,-1]) - tmax_exp_prism
tmin_anom <- as.matrix(tmin[,-1]) - tmin_exp_prism

##Plots everything to make sure it works.
par(mfrow=c(3,3))
for(i in 1:ncol(tavg_anom)){
  try(plot(density(tavg_anom[,i],na.rm=TRUE),main=colnames(tavg_anom)[i],xlim=c(-15,15)))
  try(points(density(tmax_anom[,i],na.rm=TRUE),col=2,type="l"))
  try(points(density(tmin_anom[,i],na.rm=TRUE),col=4,type="l"))
}

par(mfrow=c(1,1))
plot(tmax_anom[,177],tmax_anom[,178],pch=".")
abline(0,1,lty=2)

##Writes the PRISM seasonal averages to disk.
write.csv(days_prism,"prism_seas_2009_2015.csv",row.names=FALSE)
write.csv(days_prism_exp,"prism_seas_avg_interp_2009_2015.csv",row.names=FALSE)

##Writes raw anomalies to disk.
tavg_anom_df <- data.frame(date=prism_exp_tavg$date,tavg_anom)
tmin_anom_df <- data.frame(date=prism_exp_tavg$date,tmin_anom)
tmax_anom_df <- data.frame(date=prism_exp_tavg$date,tmax_anom)

write.csv(tavg_anom_df,"../cleaned/tavg_anomalies_2009_2015.csv",row.names=FALSE)
write.csv(tmin_anom_df,"../cleaned/tmin_anomalies_2009_2015.csv",row.names=FALSE)
write.csv(tmax_anom_df,"../cleaned/tmax_anomalies_2009_2015.csv",row.names=FALSE)


