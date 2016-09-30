##Script to merge data from the Alpine Lakes, NPS, and MesoWest datasets.
library(xts)
source("~/code/MORA_microclimate/air_temp_functions.R")
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/compiled_NPS_MesoWest/")

##Reads in hourly data.
mw_hourly <- read.csv("airtemp_MesoWest_hourly.csv")
mw_hourly$datetime <- as.POSIXct(mw_hourly$datetime,format="%Y-%m-%d %H:%M:%S",tz="Etc/GMT-8") 
mw_hourly$datetime <- as.POSIXct(format(mw_hourly$datetime,tz="America/Los_Angeles"))
mw_hourly_xts <- xts(mw_hourly[,-1],order.by=mw_hourly[,1])

al_hourly <- read.csv("airtemp_AlpineLakes_hourly.csv")
al_hourly$Datetime <- as.POSIXct(al_hourly$Datetime,format="%Y-%m-%d %H:%M:%S")
al_hourly_xts <- xts(al_hourly[,-1],order.by=al_hourly[,1])

nps_hourly <- read.csv("airtemp_NPS_hourly.csv")
nps_hourly$datetime <- as.POSIXct(nps_hourly$datetime,format="%Y-%m-%d %H:%M:%S")
nps_hourly_xts <- xts(nps_hourly[,-c(1,2)],order.by=nps_hourly$datetime)

##Prism time-series.
prism <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/PRISM_daily_1_2004_9_2015.csv") 
prism$DATE <- as.Date(prism$DATE)
prism_xts <- xts(prism[,-(1:2)],order.by=prism[,2])

##Merges MesoWest and alpine Lakes hourly data.
hourly_all <- merge(mw_hourly_xts,al_hourly_xts,nps_hourly_xts,all=TRUE)

##Computes daily statistics from hourly measurements.
data_daily_tmax_all <- apply.daily(hourly_all[,1],FUN=max)
data_daily_tmin_all <- apply.daily(hourly_all[,1],FUN=min)
data_daily_tavg_all <- apply.daily(hourly_all[,1],FUN=mean)

for(i in  2:dim(hourly_all)[[2]]){
  data_daily_tmax <- apply.daily(hourly_all[,i],FUN=max)
  data_daily_tmax_all <- merge(data_daily_tmax_all,data_daily_tmax)
  data_daily_tmin <- apply.daily(hourly_all[,i],FUN=min)
  data_daily_tmin_all <- merge(data_daily_tmin_all,data_daily_tmin)
  data_daily_tavg <- apply.daily(hourly_all[,i],FUN=mean)
  data_daily_tavg_all <- merge(data_daily_tavg_all,data_daily_tavg)
}

##Converts indices to date format.
index(data_daily_tmax_all) <- as.Date(index(data_daily_tmax_all))
index(data_daily_tmin_all) <- as.Date(index(data_daily_tmin_all))
index(data_daily_tavg_all) <- as.Date(index(data_daily_tavg_all))

##Merges with PRISM data.
data_daily_tmax_prism <- merge(data_daily_tmax_all,prism_xts[,4])
data_daily_tmin_prism <- merge(data_daily_tmin_all,prism_xts[,3])
data_daily_tavg_prism <- merge(data_daily_tavg_all,prism_xts[,1])

##Subsets to correct date range.
date_start_tmax <- which(index(data_daily_tmax_prism)==as.Date("2009-10-01"))
date_end_tmax <- which(index(data_daily_tmax_prism)==as.Date("2015-09-30"))
data_daily_tmax_df <- data.frame(date=index(data_daily_tmax_prism)[date_start_tmax:date_end_tmax],
                                 data_daily_tmax_prism[date_start_tmax:date_end_tmax])

date_start_tmin <- which(index(data_daily_tmin_prism)==as.Date("2009-10-01"))
date_end_tmin <- which(index(data_daily_tmin_prism)==as.Date("2015-09-30"))
data_daily_tmin_df <- data.frame(date=index(data_daily_tmin_prism)[date_start_tmin:date_end_tmin],
                                 data_daily_tmin_prism[date_start_tmin:date_end_tmin])

date_start_tavg <- which(index(data_daily_tavg_prism)==as.Date("2009-10-01"))
date_end_tavg <- which(index(data_daily_tavg_prism)==as.Date("2015-09-30"))
data_daily_tavg_df <- data.frame(date=index(data_daily_tavg_prism)[date_start_tavg:date_end_tavg],
                                 data_daily_tavg_prism[date_start_tavg:date_end_tavg])

data_cols_tmin <- which(colSums(!is.na(data_daily_tmin_df[,-1]),na.rm=TRUE)>30)+1
data_cols_tmax <- which(colSums(!is.na(data_daily_tmax_df[,-1]),na.rm=TRUE)>30)+1
data_cols_tavg <- which(colSums(!is.na(data_daily_tavg_df[,-1]),na.rm=TRUE)>30)+1

data_daily_tmax_df <- data_daily_tmax_df[,c(1,data_cols_tmax)]
data_daily_tmin_df <- data_daily_tmin_df[,c(1,data_cols_tmin)]
data_daily_tavg_df <- data_daily_tavg_df[,c(1,data_cols_tavg)]

##Identifies and removes outliers with a simple model.
month_fact <- as.factor(format(data_daily_tavg_df$date,format="%m"))
for (i in 2:(dim(data_daily_tavg_df)[2]-1)){
  tavg_mod <- lm(data_daily_tavg_df[,i]~data_daily_tavg_df$TEMP*month_fact,
                 na.action="na.exclude")
  tavg_resid <-scale(residuals(tavg_mod))
  plot(density(tavg_resid,na.rm=TRUE))
  tavg_out <- which(tavg_resid > 6 | tavg_resid < -6)
  if(length(tavg_out) > 0){
    data_daily_tavg_df[tavg_out,i] <- NA
  }
}
for (i in 2:(dim(data_daily_tmin_df)[2]-1)){
  tmin_mod <- lm(data_daily_tmin_df[,i]~data_daily_tmin_df$TMIN*month_fact,
                 na.action="na.exclude")
  tmin_resid <-scale(residuals(tmin_mod))
  plot(density(tmin_resid,na.rm=TRUE))
  tmin_out <- which(tmin_resid > 6 | tavg_resid < -6)
  if(length(tmin_out) > 0){
    data_daily_tmin_df[tmin_out,i] <- NA
  }
}
for (i in 2:(dim(data_daily_tmax_df)[2]-1)){
  tmax_mod <- lm(data_daily_tmax_df[,i]~data_daily_tmax_df$TMAX*month_fact,
                 na.action="na.exclude")
  tmax_resid <-scale(residuals(tmax_mod))
  plot(density(tmax_resid,na.rm=TRUE))
  tmax_out <- which(tmax_resid > 6 | tavg_resid < -6)
  if(length(tmax_out) > 0){
    data_daily_tmax_df[tmax_out,i] <- NA
  }
}

write.csv(data_daily_tmax_df,"MW_NPS_AlpineLakes_cleaned_daily_tmax.csv",row.names=FALSE)
write.csv(data_daily_tmin_df,"MW_NPS_AlpineLakes_cleaned_daily_tmin.csv",row.names=FALSE)
write.csv(data_daily_tavg_df,"MW_NPS_AlpineLakes_cleaned_daily_tavg.csv",row.names=FALSE)
