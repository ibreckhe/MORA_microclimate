## Script to check and clean formatted ibutton and HOBO air temperature data
## Author: Ian Breckheimer
## Date: 30 September 2014

## Sets up Workspace
library(xts)
library(zoo)
library(plyr)
library(data.table)
library(dplyr)
library(psych)
library(raster)
library(rgdal)

## Sets working directory and appropriate time zone.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/compiled_airtemp")
Sys.setenv(TZ="Etc/GMT-7")

#### Creates time-series figures for each file to check manually. ####

airtemp_files <- list.files(".",pattern=".csv$")

for (i in 1:length(airtemp_files)){
  
  # Prints progress to console.
  flush.console()
  print(paste("Now producing figure for ",airtemp_files[i],". File ",i," of ",
              length(airtemp_files)))
  #Reads data and converts to date and time-series format.
  data <- read.csv(airtemp_files[i],header=TRUE)
  data <- data[complete.cases(data),]
  data$datestring <- paste(data$YEAR,data$MONTH,data$DAY,sep="-")
  data$timestring <- paste(data$HOUR,data$MIN,"00",sep=":")
  data$datetime <- strptime(paste(data$datestring,data$timestring),
                            format="%Y-%m-%d %H:%M:%S")
  tempts <- xts(data$TEMP,order.by=data$datetime)
  
  # Save figure as pdf
  figdir <- "~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/figs/air/"
  plotname <- strsplit(airtemp_files[i], ".csv")[[1]]
  figpath <- paste(figdir,plotname,".pdf", sep = "")
  pdf(file = figpath, width = 10, height = 7)
  
  # Left Y axis
  par(mar = c(4, 6, 3, 6))
  plot(data$datetime,data$TEMP,type='l', main = plotname, xlab='Date', 
       ylab=expression(paste("Temperature", degree, "C")),col='blue',
       axes = T,ylim=c(-10,35),pch = 20)  
  dev.off()
}

#### Aligns all of the data to check coherence and check stats ####

## Creates a list of time-series from the files.
series_list <- list()
daytemp_fails <- c()
daydiffs <- c()
for (i in 1:length(airtemp_files)){

  flush.console()
  print(paste("Now processing",airtemp_files[i],". File ",i," of ",
              length(airtemp_files)))
  
  #Reads data and converts to date and time-series format.
  data <- read.csv(airtemp_files[i],header=TRUE)
  data <- data[complete.cases(data),]
  data$datestring <- paste(data$YEAR,data$MONTH,data$DAY,sep="-")
  data$timestring <- paste(data$HOUR,data$MIN,"00",sep=":")
  data$datetime <- as.POSIXct(paste(data$datestring,data$timestring),
                            format="%Y-%m-%d %H:%M:%S",tz="Etc/GMT-7")
  tempts <- xts(data$TEMP,order.by=data$datetime)
  series_list[[i]] <- tempts
  
  ## Checks that the day is warmer than the night, on average.
  day_data <- data[data$HOUR>=12,]
  night_data <- data[data$HOUR<12,]
  day_temp <- mean(day_data$TEMP,na.rm=TRUE)
  night_temp <- mean(night_data$TEMP,na.rm=TRUE)
  daydiff <- day_temp - night_temp
  daydiffs <- c(daydiffs,daydiff)
  daytemp_test <- day_temp <= night_temp
  daytemp_fails <- c(daytemp_fails,daytemp_test)
}
names(series_list) <- airtemp_files
sum(daytemp_fails)
cbind(airtemp_files,daydiffs)
airtemp_files[daytemp_fails]

## Creates a separate list for the daytemp fails.
series_fail <- series_list[daytemp_fails]
series_pass <- series_list[daytemp_fails==FALSE]

## Creates the plot
xlimits <- as.POSIXct(c("2010-8-15 00:00:00","2010-8-20 00:00:00"),tz="Etc/GMT-7")
plot(series_list[[1]],
     ylim=c(-20,40),
     ylab=expression(paste("Temperature, ", degree, "C")),
     xlim=xlimits,
     main='All Datasets',
     major.ticks=FALSE,
     minor.ticks=TRUE,
     auto.grid=FALSE,
     axes=TRUE)
for(i in 1:length(series_list)){
  if (i %in% c(1:26)){
    lines(y=series_list[[i]],x=index(series_list[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="black")
  }else if(i %in% c(27:33)){
    lines(y=series_list[[i]],x=index(series_list[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="yellow")
  }else if(i %in% c(34:59)){
    lines(y=series_list[[i]],x=index(series_list[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="purple")
  }else if(i %in% c(60:315)){
      lines(y=series_list[[i]],x=index(series_list[[i]]),
           ylim=c(-20,40),
           xlim=xlimits,
           main="",
           col="red")
  }else if (i %in% c(316:384)){
    lines(y=series_list[[i]],x=index(series_list[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="blue")
  }else if (i %in% c(385:397)){
    lines(y=series_list[[i]],x=index(series_list[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="green")
  }
}
axis.POSIXct(1, x=index(series_list[[i]]),labels = TRUE)

## Computes hour of the day with the highest temperature for each day.
time_max_series <- list() 

for (i in 1:length(series_list)){
  series <- series_list[[i]]
  time_max_series[[i]] <- do.call(rbind, lapply(split(series,"days"), function(x) x[which.max(x)]))
}
names(time_max_series) <- airtemp_files

##Computes circular mean hour of max temperature.
means <- c()
for (i in 1:length(time_max_series)){
  index <- index(time_max_series[[i]])
  hrs <- as.numeric(format(index,format="%H"))
  mean <- circadian.mean(hrs,hours=TRUE)
  means <- c(means,mean)
}
plot(density(means),rug=TRUE)

## Creates a logical vector indicating if the series is likely in timezone GMT
gmt <- means > 20

## colors series by inferred time-zone.
## Creates the plot
xlimits <- as.POSIXct(c("2015-4-17 00:00:00","2015-4-24 00:00:00"),tz="Etc/GMT-7")
plot(series_list[[1]],
     ylim=c(-20,40),
     ylab=expression(paste("Temperature, ", degree, "C")),
     xlim=xlimits,
     main='All Datasets',
     major.ticks=FALSE,
     minor.ticks=TRUE,
     auto.grid=FALSE,
     axes=TRUE)
for(i in 1:length(series_list)){
  if (gmt[i]){
    lines(y=series_list[[i]],x=index(series_list[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="blue")
  }else{
    lines(y=series_list[[i]],x=index(series_list[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="red")
  }
}

##Creates new time-series with corrected time zone.
series_list_tz <- list()
for (i in 1:length(airtemp_files)){
  
  flush.console()
  print(paste("Now processing",airtemp_files[i],". File ",i," of ",
              length(airtemp_files)))
  
  #Reads data and converts to date and time-series format.
  data <- read.csv(airtemp_files[i],header=TRUE)
  data <- data[complete.cases(data),]
  data$datestring <- paste(data$YEAR,data$MONTH,data$DAY,sep="-")
  data$timestring <- paste(data$HOUR,data$MIN,"00",sep=":")
  
  # inferred time-zone based on hour of max temperature.
  if(gmt[i]){
      data$datetime <- as.POSIXct(paste(data$datestring,data$timestring),
                              format="%Y-%m-%d %H:%M:%S",tz="Etc/GMT-7")
      tempts <- xts(data$TEMP,order.by=data$datetime)
  }else{
      data$datetime <- as.POSIXct(paste(data$datestring,data$timestring),
                                format="%Y-%m-%d %H:%M:%S",tz="GMT")
      tempts <- xts(data$TEMP,order.by=data$datetime)
      indexTZ(tempts) <- "Etc/GMT-7"
  }
  indexTZ(tempts) <- "GMT"
  series_list_tz[[i]] <- tempts
}

## Checks to see if it worked.
xlimits <- as.POSIXct(c("2014-8-17 00:00:00","2014-8-24 00:00:00"))
plot(series_list[[1]],
     ylim=c(-20,40),
     ylab=expression(paste("Temperature, ", degree, "C")),
     xlim=xlimits,
     main='All Datasets',
     major.ticks=FALSE,
     minor.ticks=TRUE,
     auto.grid=FALSE,
     axes=TRUE)
for(i in 1:length(series_list_tz)){
  if (gmt[i]){
    lines(y=series_list_tz[[i]],x=index(series_list_tz[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="blue")
  }else{
    lines(y=series_list_tz[[i]],x=index(series_list_tz[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col="red")
  }
}

## Makes a nice figure.
xlimits <- as.POSIXct(c("2014-10-1 00:00:00","2014-12-1 00:00:00"))
plot(series_list[[1]],
     ylim=c(-20,40),
     ylab=expression(paste("Temperature, ", degree, "C")),
     xlim=xlimits,
     main='All Datasets',
     major.ticks=FALSE,
     minor.ticks=FALSE,
     auto.grid=FALSE,
     axes=FALSE)
for(i in 1:length(series_list_tz)){
    lines(y=series_list_tz[[i]],x=index(series_list_tz[[i]]),
          ylim=c(-20,40),
          xlim=xlimits,
          main="",
          col=i,
          lwd=0.1)
}
axis.POSIXct(1, x=index(series_list[[i]]),labels = TRUE)
axis(2, at=c(-20,-10,0,10,20,30,40),labels = TRUE)

## Exports cleaned time-series
setwd("../cleaned_airtemp2")
for (i in 1:length(airtemp_files)){
  series <- series_list_tz[[i]]
  outdf <- data.frame(DATE=index(series),TEMP=series)
  write.csv(outdf,file=airtemp_files[i],row.names=FALSE)
}

## Creates metadata for clean time-series.
nfiles <- length(airtemp_files)
metadata <- data.frame(FILE=rep("NA",nfiles),
                       DATE_MIN=rep(as.POSIXct("0001-01-01"),nfiles),
                       DATE_MAX=rep(as.POSIXct("0001-01-01"),nfiles),
                       TEMP_MIN=rep(NA,nfiles),
                       TEMP_MAX=rep(NA,nfiles),
                       LOG_INT=rep(NA,nfiles))
metadata$FILE <- as.character(metadata$FILE)
for (i in 1:length(airtemp_files)){
  series <- series_list_tz[[i]]
  
  ## Gets file name
  metadata[i,1] <- airtemp_files[i]
  
  ## Gets date range.
  metadata[i,2] <- min(index(series))
  metadata[i,3] <- max(index(series))
  
  ## Gets the temp range
  metadata[i,4] <- min(series)
  metadata[i,5] <- max(series)
  
  ## Gets the logging interval.
  metadata[i,6] <- log_int <- as.numeric(index(series[2]) - index(series[1]))
}
write.csv(metadata,file="metadata3.txt",row.names=FALSE)

## Reads in cleaned metadata and adds geographic coordinates.
setwd("../cleaned")
meta_cleaned <- read.table("metadata_cleaned_10_2015.txt",header=TRUE,sep="\t")
meta_cleaned <- unique(meta_cleaned)
meta_cleaned$DATE_MIN <- as.POSIXct(as.character(meta_cleaned$DATE_MIN),format="%m/%d/%y %H:%M")
meta_cleaned$DATE_MAX <- as.POSIXct(as.character(meta_cleaned$DATE_MAX),format="%m/%d/%y %H:%M")
sensor_locs <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/sensor_locations_updated_10-12-2015.csv")
sensor_locs$X <- as.numeric(as.character(sensor_locs$X))
sensor_locs$Y <- as.numeric(as.character(sensor_locs$Y))
sensor_locs$Long <- as.numeric(as.character(sensor_locs$Long))
sensor_locs$Lat <- as.numeric(as.character(sensor_locs$Lat))
meta_cleaned$combined_1 <- paste(meta_cleaned$SITE,meta_cleaned$Sensor,sep="-")
meta_cleaned_geo <- merge(meta_cleaned,sensor_locs,by="combined_1",all.x=TRUE)

## Appends elevation and canopy cover data to metadata.
meta_cleaned_geo$Long <- as.numeric(as.character(meta_cleaned_geo$Long))
meta_cleaned_geo$Lat <- as.numeric(as.character(meta_cleaned_geo$Lat))
coordinates(meta_cleaned_geo) <-  ~Long+Lat
proj4string(meta_cleaned_geo) <- "+proj=longlat +datum=WGS84 +no_defs"
meta_cleaned_utm <- spTransform(meta_cleaned_geo,
                                CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
elev <- raster("/Volumes/ib_working/GIS/MORA_elev_3m.tif")
elev243 <- raster("/Volumes/ib_working/GIS/MORA_elev__focal_243m.tif")
canopy <- raster("/Volumes/ib_working/GIS/MORA_can_pct_focal81m.tif")
cair <- raster("/Volumes/ib_working/GIS/MORA_coldair_3m.tif")
mora_rasters <- stack(elev,elev243,canopy,cair)
meta_cleaned_elev <- extract(x=mora_rasters,y=meta_cleaned_utm,method='simple',sp=TRUE)
meta_cleaned_elev$relev243 <- meta_cleaned_elev$MORA_elev_3m - meta_cleaned_elev$MORA_elev__focal_243m
meta_cleaned_elev$subalpine <- meta_cleaned_elev$MORA_elev_3m>1400 & meta_cleaned_elev$MORA_can_pct_focal81m < 0.75
meta_cleaned_elev$forest_high <- meta_cleaned_elev$MORA_elev_3m>1400 & meta_cleaned_elev$MORA_can_pct_focal81m >= 0.75
meta_cleaned_elev$forest_low <- meta_cleaned_elev$MORA_elev_3m<=1400 & meta_cleaned_elev$MORA_can_pct_focal81m >= 0.75
meta_cleaned_elev$clearing <- meta_cleaned_elev$MORA_can_pct_focal81m <= 0.75 & meta_cleaned_elev$MORA_elev_3m<=1400
meta_cleaned_elev$ridge <- meta_cleaned_elev$MORA_coldair_3m <= mean(meta_cleaned_elev$MORA_coldair_3m,na.rm=TRUE)
meta_cleaned_elev$cove <- meta_cleaned_elev$MORA_coldair_3m > mean(meta_cleaned_elev$MORA_coldair_3m,na.rm=TRUE)


#Gets the metadata in the same order as the files.
meta_cleaned_all <- meta_cleaned_elev[order(meta_cleaned_elev$FILE),]

## Makes a quick plot of the site categories.
par(mar=c(4,4,1,1),mfrow=c(1,1))
plot(MORA_elev_3m~MORA_can_pct_focal81m,data=meta_cleaned_elev,type="n",xlab="Canopy Cover (%)",ylab="Elevation (m)")
points(jitter(meta_cleaned_elev$MORA_can_pct_focal81m[meta_cleaned_elev$subalpine],amount=0.01),
       jitter(meta_cleaned_elev$MORA_elev_3m[meta_cleaned_elev$subalpine],amount=20),col=1)
points(jitter(meta_cleaned_elev$MORA_can_pct_focal81m[meta_cleaned_elev$forest_high],amount=0.01),
       jitter(meta_cleaned_elev$MORA_elev_3m[meta_cleaned_elev$forest_high],amount=20),col=2)
points(jitter(meta_cleaned_elev$MORA_can_pct_focal81m[meta_cleaned_elev$forest_low],amount=0.01),
       jitter(meta_cleaned_elev$MORA_elev_3m[meta_cleaned_elev$forest_low],amount=20),col=3)
points(jitter(meta_cleaned_elev$MORA_can_pct_focal81m[meta_cleaned_elev$clearing],amount=0.01),
       jitter(meta_cleaned_elev$MORA_elev_3m[meta_cleaned_elev$clearing],amount=20),col=4)
legend('bottomleft',bty="n",legend=c('Subalpine','High Forest','Low Forest','Clearing/Gap'),col=c(1:4),pch=c(21))

##Plot of ridge vs. cove sites.
plot(MORA_elev_3m~MORA_can_pct_focal81m,data=meta_cleaned_elev,type="n",xlab="Canopy Cover (%)",ylab="Elevation (m)")
points(jitter(meta_cleaned_elev$MORA_can_pct_focal81m[meta_cleaned_elev$ridge],amount=0.01),
       jitter(meta_cleaned_elev$MORA_elev_3m[meta_cleaned_elev$ridge],amount=20),col=1)
points(jitter(meta_cleaned_elev$MORA_can_pct_focal81m[meta_cleaned_elev$cove],amount=0.01),
       jitter(meta_cleaned_elev$MORA_elev_3m[meta_cleaned_elev$cove],amount=20),col=2)
legend('bottomleft',bty="n",legend=c('Ridge','Cove'),col=c(1:2),pch=c(21))


## Reads the hourly data back in and associates it with the metadata.
setwd("../cleaned_airtemp")
hourly_files <- list.files(".",pattern=".csv$")
hourly_series <- list()

##Associates metadata with the file name.
files_meta <- merge(data.frame(FILE=hourly_files),meta_cleaned_all@data,all.x=TRUE)

for (i in 1:length(hourly_files)){
  print(paste("Now processing ",hourly_files[i],". File ",i ," of",length(hourly_files)))
  hourly_meta <- as.list(files_meta[i,c(1:12,18:28)])
  temp_series <- read.csv(hourly_files[i])
  dates <- as.POSIXlt(as.character(temp_series$DATE),format="%Y-%m-%d %H:%M:%S",TZ="Etc/GMT-7")
  temp_xts <- xts(temp_series$TEMP,order.by=dates,tzone="Etc/GMT-7")
  temp_xts_no_na <- temp_xts[!is.na(index(temp_xts))]
  hourly_meta$temp <- temp_xts_no_na
  names(hourly_meta) <- c("file","code","study","site","sensor","date_min","date_max",
                          "temp_min","temp_max","log_interval_hrs","long","lat","elev",
                          "elev243","relev243","cancov81","coldair_index","subalpine","forest_high",
                          "forest_low","clearing","ridge","cove","data")
  hourly_series[[i]] <- hourly_meta
}
names(hourly_series) <- hourly_files
save(hourly_series,file="hourly_series.Rdata")
