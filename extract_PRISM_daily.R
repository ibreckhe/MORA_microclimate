## Script to extract daily time-series from PRISM 4km data.
## Author: Ian Breckheimer.

## Sets up workspace.
library(raster)

## Creates a bounding box around Mt. Rainier.
ran_ext <- extent(matrix(c(-122.172553,-120.914619,46.487004,47.212064),
                         ncol=2,byrow=TRUE))

## Imports precip data as a series of raster stacks, crops the data and writes it to disk.
setwd("~/Desktop/PRISM_daily4km")
dirs <- list.dirs(".")
for (i in 2:length(dirs)){
  flush.console()
  setwd(paste("~/Desktop/PRISM_daily4km/",dirs[i],sep=""))
  files <- list.files(".",pattern=".bil$")
  for (j in 1:length(files)){
    setwd(paste("~/Desktop/PRISM_daily4km/",dirs[i],sep=""))
    print(paste("Now processing ",files[j]))
    yeardat <- raster(files[j])
    setwd("~/Desktop/PRISM_prcp_cropped")
    years <- c(".","2004","2005","2006","2007","2008","2009",
               "2010","2011","2012","2013","2014")
    filename <- files[j]
    yearcrop <- crop(yeardat,ran_ext,filename=filename,
                     datatype='INT4S', overwrite=TRUE,progress="text")
  }
}

## Imports tavg data as a series of raster stacks, crops the data and writes it to disk.
setwd("~/Desktop/PRISM_temp_4km/")
dirs <- list.dirs(".")
for (i in 2:length(dirs)){
  flush.console()
  setwd(paste("~/Desktop/PRISM_temp_4km/",dirs[i],sep=""))
  files <- list.files(".",pattern=".bil$")
  for (j in 1:length(files)){
    setwd(paste("~/Desktop/PRISM_temp_4km/",dirs[i],sep=""))
    print(paste("Now processing ",files[j]))
    yeardat <- raster(files[j])
    setwd("~/Desktop/PRISM_temp_cropped")
    years <- c(".","2004","2005","2006","2007","2008","2009",
               "2010","2011","2012","2013","2014")
    filename <- files[j]
    yearcrop <- crop(yeardat,ran_ext,filename=filename,
                     datatype='INT4S', overwrite=TRUE,progress="text")
  }
}

## Computes an average daily temperature and precip for that day across the whole study extent.
setwd("~/Desktop/PRISM_prcp_cropped/")
prcp_files <- list.files(".",pattern=".bil$")
prcp_daily_avgs <- c()
for (i in 1:length(prcp_files)){
  flush.console()
  print(paste("Now processing ",prcp_files[i]))
  rast <- raster(prcp_files[i])
  avg <- cellStats(rast,'mean')
  prcp_daily_avgs <- c(prcp_daily_avgs,avg)
}

setwd("~/Desktop/PRISM_temp_cropped/")
temp_files <- list.files(".",pattern=".bil$")
temp_daily_avgs <- c()
for (i in 1:length(temp_files)){
  flush.console()
  print(paste("Now processing ",temp_files[i]))
  rast <- raster(temp_files[i])
  avg <- cellStats(rast,'mean')
  temp_daily_avgs <- c(temp_daily_avgs,avg)
}

## Writes to disk.
days <- seq(as.Date("2004-01-01"),as.Date("2014-3-31"),by="day")
PRISM_daily <- data.frame(DATE=days,TEMP=temp_daily_avgs,PRCP=prcp_daily_avgs)
write.csv(PRISM_daily,"~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/PRISM_daily_2004_2014.csv")

## Plots the temp data.
PRISM_daily$TEMP <- xts(PRISM_daily$TEMP,order.by=PRISM_daily$DATE)
plot(PRISM_daily$TEMP,xlim=as.POSIXct(c("2007-10-1","2014-10-1")))

