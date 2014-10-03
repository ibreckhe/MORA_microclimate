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

## Sets working directory.
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/compiled_airtemp")


#### Creates time-series figures for each file to check manually. ####

airtemp_files <- list.files(".")

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
  figdir <- "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/figs/air/"
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
xlimits <- as.POSIXct(c("2009-8-17 00:00:00","2009-8-24 00:00:00"),tz="Etc/GMT-7")
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
## Creates the plot
xlimits <- as.POSIXct(c("2010-11-17 00:00:00","2010-11-24 00:00:00"))
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
## Checks to see if it worked.
## Creates the plot
xlimits <- as.POSIXct(c("2006-10-1 00:00:00","2014-10-1 00:00:00"))
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
setwd("../cleaned_airtemp")
for (i in 1:length(airtemp_files)){
  series <- series_list_tz[[i]]
  outdf <- data.frame(DATE=index(series),TEMP=series)
  write.csv(outdf,file=airtemp_files[i],row.names=FALSE)
}

## Aggregates data to daily and writes to disk.
setwd("../cleaned_dailyair")
daily_list <- list()
for (i in 1:length(airtemp_files)){
  series <- series_list_tz[[i]]
  dailymean <- apply.daily(series,FUN=function(x) mean(x))
  dailymax <- apply.daily(series,FUN=function(x) max(x))
  dailymin <- apply.daily(series,FUN=function(x) min(x))
  daily_series <- cbind(dailymean,dailymax,dailymin)
  daily_list[[i]] <- daily_series
  colnames(daily_series) <- c("TAVG","TMAX","TMIN")
  indexFormat(daily_series) <- "%Y-%m-%d"
  index <- as.character(index(daily_series),format="%Y-%m-%d")
  daily_out <- data.frame(DATE=index,daily_series)
  write.csv(daily_out,file=airtemp_files[i],row.names=FALSE)
}

## Plots daily data.
## Makes a nice figure.
## Checks to see if it worked.
## Creates the plot
xlimits <- as.POSIXct(c("2006-10-1 00:00:00","2014-10-1 00:00:00"))
plot(daily_list[[1]][,1],
     ylim=c(-20,40),
     ylab=expression(paste("Temperature, ", degree, "C")),
     xlim=xlimits,
     main='All Datasets',
     major.ticks=FALSE,
     minor.ticks=FALSE,
     auto.grid=FALSE,
     axes=FALSE)
for(i in 1:length(series_list_tz)){
  lines(y=daily_list[[i]][,1],x=index(daily_list[[i]]),
        ylim=c(-20,40),
        xlim=xlimits,
        main="",
        col=i,
        lwd=0.1)
}
axis.POSIXct(1, x=index(daily_list[[i]]),labels = TRUE)
axis(2, at=c(-20,-10,0,10,20,30,40),labels = TRUE)

