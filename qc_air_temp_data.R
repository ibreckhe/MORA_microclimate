## Script to qc and explore formatted and cleaned daily ibutton and HOBO air temperature data
## Author: Ian Breckheimer
## Date: 8 October 2014

#### Sets up Workspace and loads data ####
library(raster)
library(rgdal)
library(xts)


Sys.setenv(TZ="Etc/GMT-7")

## Reads in hourly and daily temp data
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/cleaned_airtemp/")
load("hourly_series.Rdata")

## Loads in PRISM and topoWX time-series.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
prism <- read.csv("PRISM_daily_2004_2015.csv")
topowx <- read.csv("topowx_daily_2004_2012.csv")

## Interpolate hourly series to fill missing hours.
interp_series <- hourly_series
for (i in 1:length(hourly_series)){
  print(paste("Interpolating ",hourly_series[[i]]$file,". Dataset ",i ," of",length(hourly_series)))
  temp <- hourly_series[[i]]$data
  index(temp) <- as.POSIXct(index(temp),tz="Etc/GMT-7")
  if(as.numeric(strftime(index(temp)[1],format="%M"))==0){
    min_date <- as.POSIXct(min(index(temp)),tz="Etc/GMT-7")
  }else{
    min_date <- as.POSIXct(as.Date(min(index(temp))),tz="Etc/GMT-7")
  }
  max_date <- max(index(temp))
  hours <- seq(from=min_date,
               to=max_date,by="hour")
  if(length(temp)==length(hours) & as.numeric(strftime(index(temp)[1],format="%M"))==0){
    print("Data already hourly, skipping interpolation.")
    interp_series[[i]]$temp <- temp
    interp_series[[i]]$data <- data.frame(TEMP=temp,
                                          HOURS=as.numeric(strftime(index(temp),format="%H")),
                                          MONTHS=as.numeric(strftime(index(temp),format="%m")))
  }else{
    temp_na <- cbind(temp,hours)
    temp_interp <- na.approx(object=temp_na,maxgap=4)
    interp_series[[i]]$temp <- temp_interp
    interp_series[[i]]$data <- data.frame(TEMP=temp_interp,
                                          HOURS=as.numeric(strftime(index(temp_interp),format="%H")),
                                          MONTHS=as.numeric(strftime(index(temp_interp),format="%m")))
  }
}

## Checks to make sure all the series have data.
n_measures <- rep(NA,length(hourly_series))
for (i in 1:length(hourly_series)){
  data <- hourly_series[[i]]$data
  n_measures[i] <- sum(!is.na(data))
  hourly_series[[i]]$n_measurements <- n_measures[i]
}

#### Quantify microclimate in different environments ####

## Filters out suspect sensors.
interp_series <- Filter(function(x) x$site!="TO11-STR",interp_series)

## Subsets interpolated series by biome.
subalpine <- interp_series[unlist(sapply(interp_series,"[","subalpine"))]
forest_high <- interp_series[unlist(sapply(interp_series,"[","forest_high"))]
forest_low <- interp_series[unlist(sapply(interp_series,"[","forest_low"))]
clearing <- interp_series[unlist(sapply(interp_series,"[","clearing"))]
ridge <- interp_series[unlist(sapply(interp_series,"[","ridge"))]
cove <- interp_series[unlist(sapply(interp_series,"[","cove"))]


tod.avg.fun <- function(series_list,months,hours){
  output <- matrix(NA,ncol=length(hours),nrow=length(series_list))
  for (i in 1:length(series_list)){
    if (is.null(series_list[[i]]) == FALSE){
      for (j in hours){
        sub_months <- series_list[[i]]$data[series_list[[i]]$data$MONTHS %in% months,]
        sub_hrs <- sub_months[sub_months$HOURS==j,1]
        mean <- mean(sub_hrs)
        output[i,(j+1)] <- mean
      }
    }
  }
  colnames(output) <- paste("hr_",hours,sep="")
  rownames(output) <- names(series_list)
  return(output)
}

subalpine_summer <- tod.avg.fun(subalpine,months=c(6:8),hours=c(0:23))
forest_high_summer <- tod.avg.fun(forest_high,months=c(6:8),hours=c(0:23))
forest_low_summer <- tod.avg.fun(forest_low,months=c(6:8),hours=c(0:23))
clearing_summer <- tod.avg.fun(clearing,months=c(6:8),hours=c(0:23))
ridge_summer <- tod.avg.fun(ridge,months=c(6:8),hours=c(0:23))
cove_summer <- tod.avg.fun(cove,months=c(6:8),hours=c(0:23))


subalpine_fall <- tod.avg.fun(subalpine,months=c(9:11),hours=c(0:23))
forest_high_fall <- tod.avg.fun(forest_high,months=c(9:11),hours=c(0:23))
forest_low_fall <- tod.avg.fun(forest_low,months=c(9:11),hours=c(0:23))
clearing_fall <- tod.avg.fun(clearing,months=c(9:11),hours=c(0:23))
ridge_fall <- tod.avg.fun(ridge,months=c(9:11),hours=c(0:23))
cove_fall <- tod.avg.fun(cove,months=c(9:11),hours=c(0:23))


subalpine_winter <- tod.avg.fun(subalpine,months=c(12,1,2),hours=c(0:23))
forest_high_winter <- tod.avg.fun(forest_high,months=c(12,1,2),hours=c(0:23))
forest_low_winter <- tod.avg.fun(forest_low,months=c(12,1,2),hours=c(0:23))
clearing_winter <- tod.avg.fun(clearing,months=c(12,1,2),hours=c(0:23))
ridge_winter <- tod.avg.fun(ridge,months=c(12,1,2),hours=c(0:23))
cove_winter <- tod.avg.fun(cove,months=c(12,1,2),hours=c(0:23))


subalpine_spring <- tod.avg.fun(subalpine,months=c(3:5),hours=c(0:23))
forest_high_spring <- tod.avg.fun(forest_high,months=c(3:5),hours=c(0:23))
forest_low_spring<- tod.avg.fun(forest_low,months=c(3:5),hours=c(0:23))
clearing_spring <- tod.avg.fun(clearing,months=c(3:5),hours=c(0:23))
ridge_spring <- tod.avg.fun(ridge,months=c(3:5),hours=c(0:23))
cove_spring <- tod.avg.fun(cove,months=c(3:5),hours=c(0:23))


## Gets colmeans.
subalpine_summer_colmeans <- colMeans(subalpine_summer,na.rm=TRUE)
forest_high_summer_colmeans <- colMeans(forest_high_summer,na.rm=TRUE)
forest_low_summer_colmeans <- colMeans(forest_low_summer,na.rm=TRUE)
clearing_summer_colmeans <- colMeans(clearing_summer,na.rm=TRUE)
ridge_summer_colmeans <- colMeans(ridge_summer,na.rm=TRUE)
cove_summer_colmeans <- colMeans(cove_summer,na.rm=TRUE)


subalpine_fall_colmeans <- colMeans(subalpine_fall,na.rm=TRUE)
forest_high_fall_colmeans <- colMeans(forest_high_fall,na.rm=TRUE)
forest_low_fall_colmeans <- colMeans(forest_low_fall,na.rm=TRUE)
clearing_fall_colmeans <- colMeans(clearing_fall,na.rm=TRUE)
ridge_fall_colmeans <- colMeans(ridge_fall,na.rm=TRUE)
cove_fall_colmeans <- colMeans(cove_fall,na.rm=TRUE)


subalpine_winter_colmeans <- colMeans(subalpine_winter,na.rm=TRUE)
forest_high_winter_colmeans <- colMeans(forest_high_winter,na.rm=TRUE)
forest_low_winter_colmeans <- colMeans(forest_low_winter,na.rm=TRUE)
clearing_winter_colmeans <- colMeans(clearing_winter,na.rm=TRUE)
ridge_winter_colmeans <- colMeans(ridge_winter,na.rm=TRUE)
cove_winter_colmeans <- colMeans(cove_winter,na.rm=TRUE)


subalpine_spring_colmeans <- colMeans(subalpine_spring,na.rm=TRUE)
forest_high_spring_colmeans <- colMeans(forest_high_spring,na.rm=TRUE)
forest_low_spring_colmeans <- colMeans(forest_low_spring,na.rm=TRUE)
clearing_spring_colmeans <- colMeans(clearing_spring,na.rm=TRUE)
ridge_spring_colmeans <- colMeans(ridge_spring,na.rm=TRUE)
cove_spring_colmeans <- colMeans(cove_spring,na.rm=TRUE)


##Plots hourly averages for each sensor.
par(mfrow=c(1,1),mar=c(2,2,2,2),oma=c(1,1,1,1))
plot(0:23,subalpine_summer_colmeans,type="n", ylim=c(-5,25),axes=TRUE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
for (i in 1:dim(forest_low_summer)[2]){
  lines(0:23,forest_low_summer[i,],col=1)
}
for (i in 1:dim(clearing_summer)[2]){
  lines(0:23,clearing_summer[i,],col=2)
}
for (i in 1:dim(forest_high_summer)[2]){
  lines(0:23,forest_high_summer[i,],col=3)
}
for (i in 1:dim(subalpine_summer)[2]){
  lines(0:23,subalpine_summer[i,],col=4)
}

##Makes a plot of average hourly conditions by season.
pdf("../results/microclimate_biome_seasonal2.pdf",width=5,height=5)
par(mfrow=c(2,2),mar=c(1,1,1,1),oma=c(4,4,0,0))
plot(0:23,subalpine_winter_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25))
axis(1,at=c(seq(2,22,by=4)),tick=TRUE,labels=FALSE)
abline(h=0,col="blue",lwd=0.5)
lines(0:23,forest_low_winter_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,clearing_winter_colmeans,type="l",col="black",lwd=3,lty=4)
lines(0:23,forest_high_winter_colmeans,type="l",col="grey60",lwd=2,lty=1)
lines(0:23,subalpine_winter_colmeans,type="l",col="grey60",lwd=3,lty=2)
text(x=0,y=22,pos=4,labels="Winter (DJF)",font=2,cex=1.1)
legend(2,20,bty="n",legend=c("Low-elev. Forest","Clearing","High-elev. Forest","Subalpine"),lty=c(1,4,1,2),
       lwd=c(2,3,2,3),col=c("black","black","grey60","grey60"),y.intersp=0.8)

plot(0:23,subalpine_spring_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25),tick=TRUE,labels=FALSE)
axis(1,at=c(seq(2,22,by=4)),tick=TRUE,labels=FALSE)
abline(h=0,col="blue",lwd=0.5)
lines(0:23,forest_low_spring_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,clearing_spring_colmeans,type="l",col="black",lwd=3,lty=4)
lines(0:23,forest_high_spring_colmeans,type="l",col="grey60",lwd=2,lty=1)
lines(0:23,subalpine_spring_colmeans,type="l",col="grey60",lwd=3,lty=2)
text(x=0,y=22,pos=4,labels="Spring (MAM)",font=2,cex=1.1)

plot(0:23,subalpine_winter_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25),labels=c(-5,0,5,10,15,20,25))
axis(1,at=c(seq(2,22,by=4)))
abline(h=0,col="blue",lwd=0.5)
lines(0:23,forest_low_summer_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,clearing_summer_colmeans,type="l",col="black",lwd=3,lty=4)
lines(0:23,forest_high_summer_colmeans,type="l",col="grey60",lwd=2,lty=1)
lines(0:23,subalpine_summer_colmeans,type="l",col="grey60",lwd=3,lty=2)
text(x=0,y=22,pos=4,labels="Summer (JJA)",font=2,cex=1.1)

plot(0:23,subalpine_summer_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25),tick=TRUE,labels=FALSE)
axis(1,at=c(seq(2,22,by=4)))
abline(h=0,col="blue",lwd=0.5)
lines(0:23,forest_low_fall_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,clearing_fall_colmeans,type="l",col="black",lwd=3,lty=4)
lines(0:23,forest_high_fall_colmeans,type="l",col="grey60",lwd=2,lty=1)
lines(0:23,subalpine_fall_colmeans,type="l",col="grey60",lwd=3,lty=2)
text(x=0,y=22,pos=4,labels="Fall (SON)",font=2,cex=1.1)

mtext(text="Time of Day (Hrs.)",side=1,outer=TRUE,padj=2)
mtext(text="Air Temperature (C)",side=2,outer=TRUE,padj=-2)

dev.off()

##Makes a dataset to share with Janneke.
hourly_microclim_biome_seasonal <- data.frame(habitat=c(rep("Low-elev. Forest",4),
                                                        rep("Clearing / Gap",4),
                                                        rep("High-elev. Forest",4),
                                                        rep("Subalpine",4)),
                                              season=rep(c("Winter (DJF)","Spring (MAM)",
                                                           "Summer (JJA)","Fall (SON)"),4))

hrly_data <- rbind(forest_low_winter_colmeans,
                   forest_low_spring_colmeans,
                   forest_low_summer_colmeans,
                   forest_low_fall_colmeans,
                   clearing_winter_colmeans,
                   clearing_spring_colmeans,
                   clearing_summer_colmeans,
                   clearing_fall_colmeans,
                   forest_high_winter_colmeans,
                   forest_high_spring_colmeans,
                   forest_high_summer_colmeans,
                   forest_high_fall_colmeans,
                   subalpine_winter_colmeans,
                   subalpine_spring_colmeans,
                   subalpine_summer_colmeans,
                   subalpine_fall_colmeans)                                              

biome_season_hrly <- cbind(hourly_microclim_biome_seasonal,
                           hrly_data)     
write.csv(biome_season_hrly,"../cleaned/biome_season_hrly_summary.csv",row.names=FALSE)


##Makes a similar plot separating ridges and coves.
pdf("../results/microclimate_ridg_cove_seasonal.pdf",width=5,height=5)
par(mfrow=c(2,2),mar=c(1,1,1,1),oma=c(4,4,0,0))
plot(0:23,subalpine_winter_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25))
axis(1,at=c(seq(2,22,by=4)),tick=TRUE,labels=FALSE)
abline(h=0,col="blue",lwd=0.5)
lines(0:23,ridge_winter_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,cove_winter_colmeans,type="l",col="black",lwd=3,lty=4)
text(x=0,y=22,pos=4,labels="Winter (DJF)",font=2,cex=1.1)
legend(2,20,bty="n",legend=c("Ridge","Cove"),lty=c(1,4),
       lwd=c(2,3),col=c("black","black"),y.intersp=0.8)

plot(0:23,subalpine_spring_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25),tick=TRUE,labels=FALSE)
axis(1,at=c(seq(2,22,by=4)),tick=TRUE,labels=FALSE)
abline(h=0,col="blue",lwd=0.5)
lines(0:23,ridge_spring_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,cove_spring_colmeans,type="l",col="black",lwd=3,lty=4)
text(x=0,y=22,pos=4,labels="Spring (MAM)",font=2,cex=1.1)

plot(0:23,subalpine_winter_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25),labels=c(-5,0,5,10,15,20,25))
axis(1,at=c(seq(2,22,by=4)))
abline(h=0,col="blue",lwd=0.5)
lines(0:23,ridge_summer_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,cove_summer_colmeans,type="l",col="black",lwd=3,lty=4)
text(x=0,y=22,pos=4,labels="Summer (JJA)",font=2,cex=1.1)

plot(0:23,subalpine_summer_colmeans,type="n", ylim=c(-5,25),axes=FALSE,
     xlab="Time of Day",ylab="Air Temp. (C)",lwd=2.5,lty=3)
box()
axis(2,at=c(-5,0,5,10,15,20,25),tick=TRUE,labels=FALSE)
axis(1,at=c(seq(2,22,by=4)))
abline(h=0,col="blue",lwd=0.5)
lines(0:23,ridge_fall_colmeans,type="l",col="black",lwd=2,lty=1)
lines(0:23,cove_fall_colmeans,type="l",col="black",lwd=3,lty=4)
text(x=0,y=22,pos=4,labels="Fall (SON)",font=2,cex=1.1)

mtext(text="Time of Day (Hrs.)",side=1,outer=TRUE,padj=2)
mtext(text="Air Temperature (C)",side=2,outer=TRUE,padj=-2)

dev.off()

#### Join sensors at the same location into a single series.

## Creates vectors of sites and hours
sites <- unique(unlist(sapply(interp_series,'[[',"code")))
hours_all <- seq(as.POSIXct("2006-9-1 01:00:00",tz="Etc/GMT-7"),
                 as.POSIXct("2015-10-1 23:00:00"),by="hour")
empty_xts <- xts(rep(NA,length(hours_all)),order.by=hours_all)

## Orders data list by minimum date.
mindates <- order(as.POSIXct(unlist(sapply(interp_series,'[[',"date_min")),
                             origin="1970-01-01 00:00.00 UTC"))
hourly_ordered <- interp_series[mindates]

## Concantates all of the records from each site.
site_series <- list()
for (i in 1:length(sites)){
  series_site <- Filter(function(x) x$code==sites[i],hourly_ordered)
  series_merged <- xts(NA,order.by=as.POSIXct("1970-01-01 00:00:00"))
  for (j in 1:length(series_site)){
    series_merged <- c(series_merged,series_site[[j]]$temp)
  }
  series_merged <- series_merged[-1,]
  series_merged <- cbind(empty_xts,series_merged)[,2]
  site_series[[i]] <- series_merged
}
names(site_series) <- sites

## Plots all merged series.
setwd("../figs/air_merged")
for (i in 1:length(site_series)){
  pdf(paste(names(site_series[i]),".pdf",sep=""),width=10,height=5)
  plot(site_series[[i]])
  dev.off()
}

## Reduces data to daily.
site_series_tmax <- list()
for (i in 1:length(site_series)){
  site_series_tmax[[i]] <- apply.daily(site_series[[i]],FUN=max,na.rm=FALSE)
}
names(site_series_tmax) <- sites
site_series_tmin <- list()
for (i in 1:length(site_series)){
  site_series_tmin[[i]] <- apply.daily(site_series[[i]],FUN=min,na.rm=FALSE)
}
names(site_series_tmin) <- sites
site_series_tavg <- list()
for (i in 1:length(site_series)){
  site_series_tavg[[i]] <- apply.daily(site_series[[i]],FUN=mean,na.rm=FALSE)
}
names(site_series_tavg) <- sites

## Puts everything in a data frame.
days_all <- seq(as.POSIXct("2006-9-1 23:00:00",tz="Etc/GMT-7"),
                 as.POSIXct("2015-10-1 23:00:00"),by="day")
alldat_tmax <- data.frame(DATE=as.Date(days_all))
for (i in 1:length(site_series_tmax)){
  alldat_tmax <- cbind(alldat_tmax,data.frame(round(site_series_tmax[[i]],digits=3)))
}
colnames(alldat_tmax) <- c("DATE",sites)
alldat_tmin <- data.frame(DATE=as.Date(days_all))
for (i in 1:length(site_series_tmin)){
  alldat_tmin <- cbind(alldat_tmin,data.frame(round(site_series_tmin[[i]],digits=3)))
}
colnames(alldat_tmin) <- c("DATE",sites)
alldat_tavg <- data.frame(DATE=as.Date(days_all))
for (i in 1:length(site_series_tavg)){
  alldat_tavg <- cbind(alldat_tavg,data.frame(round(site_series_tavg[[i]],digits=3)))
}
colnames(alldat_tavg) <- c("DATE",sites)

## writes data to disk.
setwd("../")
setwd("../cleaned_dailyair/")
write.csv(alldat_tmax,"alldat_daily_tmax_2006_2015.csv",row.names=FALSE)
write.csv(alldat_tmin,"alldat_daily_tmin_2006_2015.csv",row.names=FALSE)
write.csv(alldat_tavg,"alldat_daily_tavg_2006_2015.csv",row.names=FALSE)


## Plots daily data in a ridiculously large plot.
pdf("../figs/daily_temp_ultrawide_2015.pdf",width=50,height=6)
xlimits <- as.POSIXct(c("2009-1-1 00:00:00","2015-6-30 00:00:00"))
plot(site_series_tavg[[1]],
     ylim=c(-20,30),
     ylab=expression(paste("Avg. Air Temperature, ", degree, "C")),
     xlim=xlimits,
     main='All Datasets',
     major.ticks=FALSE,
     minor.ticks=FALSE,
     auto.grid=FALSE,
     axes=FALSE,
     lwd=0.3)
for(i in 1:length(site_series_tavg)){
  lines(y=site_series_tavg[[i]],x=site_series_tavg[[i]],
        ylim=c(-20,40),
        xlim=xlimits,
        main="",
        col=i,
        lwd=0.3)
}
xlabels <- seq(as.POSIXct("2008-9-30 00:00:00"),as.POSIXct("2015-6-30 00:00:00"),by='month')
xlabels2 <- seq(as.POSIXct("2009-1-1 00:00:00"),as.POSIXct("2015-1-1 00:00:00"),by='year')
axis.POSIXct(side=1, at=xlabels,format="%b",labels = TRUE)
axis.POSIXct(side=1,at=xlabels2,format="%Y",tck=1,labels = TRUE,hadj=-9,padj=1.8)
axis(side=2, at=c(-20,-10,0,10,20,30,40),labels = TRUE)

##Adds PRISM daily averages.
PRISM_avg <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/PRISM_daily_2004_2015.csv")
PRISM_avg$DATE <- as.POSIXct(PRISM_avg$DATE)
PRISM_avg$TEMP <- xts(PRISM_avg$TEMP,order.by=PRISM_avg$DATE)
PRISM_avg <- PRISM_avg[PRISM_avg$DATE >= "2008-9-30",]
lines(PRISM_avg$TEMP,col=1,lwd=2)

dev.off()

## Gets the number of sensors for each day.
sensor_nums <- apply(alldat_tavg[,-1],FUN=function(x) sum(!is.na(x)),MARGIN=1)
mean_tmax <- apply(alldat_tmax[,-1],FUN=function(x) mean(x,na.rm=TRUE),MARGIN=1)
mean_tmin <- apply(alldat_tmin[,-1],FUN=function(x) mean(x,na.rm=TRUE),MARGIN=1)
mean_tavg <- apply(alldat_tavg[,-1],FUN=function(x) mean(x,na.rm=TRUE),MARGIN=1)

## gets metadata by site.
metadata <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/sensor_locations_updated_10-12-2015.csv")
sites <- unique(unlist(sapply(hourly_series,'[[',"code")))
metadata_merged <- merge(data.frame(sites=sites),metadata,by.x="sites",
                         by.y="combined_1",keep.x=TRUE)[,1:7]
metadata_merged$X <- as.numeric(as.character(metadata_merged$X))
metadata_merged$Y <- as.numeric(as.character(metadata_merged$Y))
coordinates(metadata_merged) <- ~X+Y
proj4string(metadata_merged) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

## transforms to UTM
metadata_utm <- spTransform(metadata_merged,
                            CRS("+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

## Brings in the raster data.
can_pct <- raster("/Volumes/ib_working/GIS/MORA_can_pct_3m.tif")
can_pct_9m <- raster("/Volumes/ib_working/GIS/MORA_can_pct_focal9m.tif")
can_pct_81m <- raster("/Volumes/ib_working/GIS/MORA_can_pct_focal81m.tif")
can_vol_81m <- raster("/Volumes/ib_working/GIS/MORA_can_vol_81m.tif")
elev <- raster("/Volumes/ib_working/GIS/MORA_elev_3m.tif")
elev_9m <- raster("/Volumes/ib_working/GIS/MORA_elev_focal9m.tif")
elev_81m <- raster("/Volumes/ib_working/GIS/MORA_elev_focal81m.tif")
elev_243m <- raster("/Volumes/ib_working/GIS/MORA_elev__focal_243m.tif")
elev_729m <- raster("/Volumes/ib_working/GIS/MORA_elev__focal_729m.tif")
srad_yearsum_9m <- raster("/Volumes/ib_working/GIS/MORA_srad_yearsum_9m.tif")
coldair_index <- raster("/Volumes/ib_working/GIS/MORA_coldair_index.tif")
dry_index_9m <- raster("/Volumes/ib_working/GIS/MORA_dry_index_9m_precip.tif")
dry_index_81m <- raster("/Volumes/ib_working/GIS/MORA_dry_index_81m_precip.tif")

rast_stack <- stack(can_pct,can_pct_9m,can_pct_81m,elev,
elev_9m,elev_81m,elev_243m,elev_729m,can_vol_81m,srad_yearsum_9m,coldair_index,
dry_index_9m,dry_index_81m)

## Samples each raster at the sensor locations.
metadata_utm$Lon <- coordinates(metadata_merged)[,1]
metadata_utm$Lat <- coordinates(metadata_merged)[,2]
metadata_utm$UTM_X<- coordinates(metadata_utm)[,1]
metadata_utm$UTM_Y <- coordinates(metadata_utm)[,2]
metadata_covars <- extract(rast_stack,metadata_utm,method="simple",fun=mean,sp=TRUE)
colnames(metadata_covars@data) <- c("location","study","site","sub_site1","sub_site2","lon","lat",
                                    "utmx","utmy","ccov","ccov_9m","ccov_81m","elev","elev_9m",
                                    "elev_81m","elev_243m","elev_729m","cvol_81m","srad_sum_9m","cold_ind_9m",
                                    "dry_ind_9m","dry_ind_81m")
metadata_covars$relev_9m <- metadata_covars$elev - metadata_covars$elev_9m
metadata_covars$relev_81m <- metadata_covars$elev - metadata_covars$elev_81m
metadata_covars$relev_243m <- metadata_covars$elev - metadata_covars$elev_243m
metadata_covars$relev_729m <- metadata_covars$elev - metadata_covars$elev_729m

## Writes the metadata to disk.
write.csv(metadata_covars@data,"airtemp_site_metadata_2015.csv",row.names=FALSE)

