##Script to clean and format mesoWest meteorological data.

library(dplyr)
library(xts)
source("~/code/MORA_microclimate/air_temp_functions.R")

setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/MesoWest/")

files <- list.files(".")
airtemp_files <- grep("air\\_temperature",files,value=TRUE)

hours <- seq(as.POSIXct("2009-01-01 00:00:00 GMT",tz="GMT"),as.POSIXct("2015-12-31 00:00:00 GMT",tz="GMT"),by="hour")
temps <- rep(NA,length(hours))
hour_ts <- xts(temps,order.by=hours,tzone="GMT")

##Loops through each file and cleans them.
for(i in 1:length(airtemp_files)){
  print(paste("Processing file" ,airtemp_files[i]))
  tdata <- read.table(airtemp_files[i],header=FALSE,sep=",",skip=8)
  colnames(tdata) <- c("Station","Date","Celcius")
  tdata$Datetime <- as.POSIXct(tdata$Date,format="%Y-%m-%dT%H:%M:%SZ",tz="GMT")
  tdata_ts <- xts(tdata$Celcius,order.by=tdata$Datetime,tzone="GMT")
  tdata_merged <- merge(hour_ts,tdata_ts,all=TRUE)
  tdata_merged[!is.na(tdata_merged[,2]),1] <- tdata_merged[!is.na(tdata_merged[,2]),2]
  tdata_all <- tdata_merged[,1]
  tdata_all_filled <- na.spline(tdata_all,maxgap=2,na.rm=FALSE)
  tdata_filled_df <- data.frame(Date=index(tdata_all_filled),Temp=tdata_all_filled)
  tdata_filled_df$hr <- as.numeric(format(tdata_filled_df$Date,format="%M"))
  tdata_filled_hrly <- tdata_filled_df[tdata_filled_df$hr==0,c(1,2)]
  tdata_filled_hrly_ts <- xts(tdata_filled_hrly$hour_ts,order.by=tdata_filled_hrly$Date)
  tdata_filled_hrly_ts2 <- remove_spikes(tdata_filled_hrly_ts,ts_lag=2,thresh=c(-5,5))
  tdata_filled_hrly_ts2[tdata_filled_hrly_ts2 > 40 | tdata_filled_hrly_ts2 < -30] <- NA
  tdata_filled_hrly_ts3 <- na.spline(tdata_filled_hrly_ts2,maxgap=2,na.rm=FALSE)
  tdata_cleaned <- data.frame(datetime=index(tdata_filled_hrly_ts3),Temp_C=as.numeric(tdata_filled_hrly_ts3))
  tdata_cleaned$datetime <- format(tdata_cleaned$datetime, tz="Etc/GMT-8",usetz=TRUE)
  write.csv(tdata_cleaned,paste("../Mesowest_cleaned/",airtemp_files[i],sep=""),row.names=FALSE)
}

##Plots all time-series to check temporal alignment.
setwd("../MesoWest_cleaned/")
cleaned_files <- list.files(".",pattern="air\\_temperature")

cdata_dates <- as.POSIXct(read.csv(cleaned_files[i])$datetime,tz="Etc/GMT-8")
empty <- rep(NA,length(cdata_dates))
cdata_all <- xts(empty,order.by=cdata_dates)

plot.xts(cdata_all,xlim=c(as.POSIXct("2011-04-30"),as.POSIXct("2011-05-28")),ylim=c(-10,20))
for (i in 1:length(cleaned_files)){
  cdata <- read.csv(cleaned_files[i])
  stationname <- strsplit(cleaned_files[i],split="\\_")[[1]][1]
  cdata_ts <- xts(cdata[,2],order.by=as.POSIXct(cdata$datetime,tz="Etc/GMT-8"))
  names(cdata_ts) <- stationname
  colnames(cdata) <- c("datetime",stationname)
  points(cdata_ts,col=i,lwd=0.5,type="l")
  cdata_all <- merge(cdata_all,cdata_ts)
}
cdata_all_df <- data.frame(datetime=index(cdata_all),cdata_all)
cdata_all_df <- cdata_all_df[,-2]
cdata_all <- cdata_all[,-1]
write.csv(cdata_all_df,"all_temp_data_hourly.csv",row.names=FALSE)

##Computes daily statistics from hourly measurements.
data_daily_tmax_all <- apply.daily(cdata_all[,1],FUN=max)
data_daily_tmin_all <- apply.daily(cdata_all[,1],FUN=min)
data_daily_tavg_all <- apply.daily(cdata_all[,1],FUN=mean)

for(i in  2:dim(cdata_all)[[2]]){
  data_daily_tmax <- apply.daily(cdata_all[,i],FUN=max)
  data_daily_tmax_all <- merge(data_daily_tmax_all,data_daily_tmax)
  data_daily_tmin <- apply.daily(cdata_all[,i],FUN=min)
  data_daily_tmin_all <- merge(data_daily_tmin_all,data_daily_tmin)
  data_daily_tavg <- apply.daily(cdata_all[,i],FUN=mean)
  data_daily_tavg_all <- merge(data_daily_tavg_all,data_daily_tavg)
}

#Compares each station with regional PRISM time-series.
prism <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/PRISM_daily_1_2004_9_2015.csv") 
prism$DATE <- as.Date(prism$DATE)
prism_tmax_ts <- xts(prism$TMAX,prism$DATE)
prism_tmin_ts <- xts(prism$TMIN,prism$DATE)
prism_tavg_ts <- xts(prism$TEMP,prism$DATE)

index(data_daily_tavg_all) <- as.Date(index(data_daily_tavg_all))
index(data_daily_tmin_all) <- as.Date(index(data_daily_tmin_all))
index(data_daily_tmax_all) <- as.Date(index(data_daily_tmax_all))
data_daily_tavg_prism <- merge(data_daily_tavg_all,prism_tavg_ts,all=c(TRUE,FALSE))
data_daily_tmin_prism <- merge(data_daily_tmin_all,prism_tavg_ts,all=c(TRUE,FALSE))
data_daily_tmax_prism <- merge(data_daily_tmax_all,prism_tavg_ts,all=c(TRUE,FALSE))

start <- which(index(data_daily_tmax_prism)==as.Date("2009-10-01"))
end <- which(index(data_daily_tmax_prism)==as.Date("2015-09-30"))

daily_tmax_df <- data.frame(datetime=as.Date(index(data_daily_tmax_prism[start:end,-49])),data_daily_tmax_prism[start:end,-49])
daily_tmin_df <- data.frame(datetime=as.Date(index(data_daily_tmin_prism[start:end,-49])),data_daily_tmin_prism[start:end,-49])
daily_tavg_df <- data.frame(datetime=as.Date(index(data_daily_tavg_prism[start:end,-49])),data_daily_tavg_prism[start:end,-49])


write.csv(daily_tmax_df,"MW_data_daily_tmax.csv",row.names=FALSE)
write.csv(daily_tmin_df,"MW_data_daily_tmin.csv",row.names=FALSE)
write.csv(daily_tmin_df,"MW_data_daily_tavg.csv",row.names=FALSE)


##Brings in site metadata.
metadata <- read.csv("../MesoWest/MW_metadata_cleanedlocs_40km.csv")
metadata$elev <- metadata$ELEVATION * 0.3048

##Estimates disparity and sensitivity for each station for each season.
start <- which(index(data_daily_tmax_prism)==as.Date("2009-10-01"))
end <- which(index(data_daily_tmax_prism)==as.Date("2015-09-30"))


indices_tmax <- boot_indices(regional_series=data_daily_tmax_prism[start:end,49],
                             local_frame=data_daily_tmax_prism[start:end,-49],
                             site_names=names(data_daily_tmax_prism)[-49],
                             nboot=1000)
colnames(indices_tmax)[-1] <- paste("tmax_",colnames(indices_tmax[-1]),sep="")

indices_tmin <- boot_indices(regional_series=data_daily_tmin_prism[start:end,49],
                             local_frame=data_daily_tmin_prism[start:end,-49],
                             site_names=names(data_daily_tmin_prism)[-49],
                             nboot=1000)
colnames(indices_tmin)[-1] <- paste("tmin_",colnames(indices_tmin[-1]),sep="")

indices_tavg <- boot_indices(regional_series=data_daily_tavg_prism[start:end,49],
                             local_frame=data_daily_tavg_prism[start:end,-49],
                             site_names=names(data_daily_tavg_prism)[-49],
                             nboot=1000)
colnames(indices_tavg)[-1] <- paste("tavg_",colnames(indices_tavg[-1]),sep="")

##Merges with metadata.
metadata_indices <- merge(metadata,indices_tavg,by.x="STID",by.y="site",all.x=TRUE)
metadata_indices <- merge(metadata_indices,indices_tmin,by.x="STID",by.y="site",all.x=TRUE)
metadata_indices <- merge(metadata_indices,indices_tmax,by.x="STID",by.y="site",all.x=TRUE)

##Visually checks daily summaries against elevation.
par(mfrow=c(2,3))
plot(tmax_disp~elev,data=metadata_indices)
abline(lm(tmax_disp~elev,data=metadata_indices),lty=2)
plot(tmin_disp~elev,data=metadata_indices)
abline(lm(tmin_disp~elev,data=metadata_indices),lty=2)
plot(tavg_disp~elev,data=metadata_indices)
abline(lm(tavg_disp~elev,data=metadata_indices),lty=2)

plot(tmax_sens~elev,data=metadata_indices)
abline(lm(tmax_sens~elev,data=metadata_indices),lty=2)
plot(tmin_sens~elev,data=metadata_indices)
abline(lm(tmin_sens~elev,data=metadata_indices),lty=2)
plot(tavg_sens~elev,data=metadata_indices)
abline(lm(tavg_sens~elev,data=metadata_indices),lty=2)

summary(lm(tmin_disp~elev+LONGITUDE,data=metadata_indices))

