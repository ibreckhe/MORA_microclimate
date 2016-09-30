## Script to extract monthly time-series from PRISM 4km data.
## Author: Ian Breckheimer.

## Sets up workspace.
library(raster)
library(xts)

## Creates a bounding box around Mt. Rainier.
ran_ext <- extent(matrix(c(-122.172553,-120.914619,46.487004,47.212064),
                         ncol=2,byrow=TRUE))

## Imports precip data as a series of raster stacks, crops the data and writes it to disk.
setwd("~/Desktop/PRISM_prcp")
dirs <- list.dirs(".")
for (i in 2:length(dirs)){
  flush.console()
  setwd(paste("~/Desktop/PRISM_prcp/",dirs[i],sep=""))
  files <- list.files(".",pattern=".bil$")
  for (j in 1:length(files)){
    setwd(paste("~/Desktop/PRISM_prcp/",dirs[i],sep=""))
    print(paste("Now processing ",files[j]))
    yeardat <- raster(files[j])
    setwd("~/Desktop/PRISM_prcp_cropped")
    filename <- files[j]
    yearcrop <- crop(yeardat,ran_ext,filename=filename,
                     datatype='INT4S', overwrite=TRUE,progress="text")
  }
}

## Imports tavg data as a series of raster stacks, crops the data and writes it to disk.
setwd("~/Desktop/PRISM_temp/")
dirs <- list.dirs(".")
for (i in 2:length(dirs)){
  flush.console()
  setwd(paste("~/Desktop/PRISM_temp/",dirs[i],sep=""))
  files <- list.files(".",pattern=".bil$")
  for (j in 1:length(files)){
    setwd(paste("~/Desktop/PRISM_temp/",dirs[i],sep=""))
    print(paste("Now processing ",files[j]))
    yeardat <- raster(files[j])
    setwd("~/Desktop/PRISM_temp_cropped")
    filename <- files[j]
    yearcrop <- crop(yeardat,ran_ext,filename=filename,
                     datatype='INT4S', overwrite=TRUE,progress="text")
  }
}

## Imports tmin data as a series of raster stacks, crops the data and writes it to disk.
setwd("~/Desktop/PRISM_tmin/")
dirs <- list.dirs(".")
for (i in 2:length(dirs)){
  flush.console()
  setwd(paste("~/Desktop/PRISM_tmin/",dirs[i],sep=""))
  files <- list.files(".",pattern=".bil$")
  for (j in 1:length(files)){
    setwd(paste("~/Desktop/PRISM_tmin/",dirs[i],sep=""))
    print(paste("Now processing ",files[j]))
    yeardat <- raster(files[j])
    setwd("~/Desktop/PRISM_tmin_cropped")
    filename <- files[j]
    yearcrop <- crop(yeardat,ran_ext,filename=filename,
                     datatype='INT4S', overwrite=TRUE,progress="text")
  }
}

## Imports tavg data as a series of raster stacks, crops the data and writes it to disk.
setwd("~/Desktop/PRISM_tmax/")
dirs <- list.dirs(".")
for (i in 2:length(dirs)){
  flush.console()
  setwd(paste("~/Desktop/PRISM_tmax/",dirs[i],sep=""))
  files <- list.files(".",pattern=".bil$")
  for (j in 1:length(files)){
    setwd(paste("~/Desktop/PRISM_tmax/",dirs[i],sep=""))
    print(paste("Now processing ",files[j]))
    yeardat <- raster(files[j])
    setwd("~/Desktop/PRISM_tmax_cropped")
    filename <- files[j]
    yearcrop <- crop(yeardat,ran_ext,filename=filename,
                     datatype='INT4S', overwrite=TRUE,progress="text")
  }
}

## Computes an average daily temperature, tmax, tmin and precip for that day across the whole study extent.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw")

setwd("./PRISM_prcp_cropped/")
prcp_files <- list.files(".",pattern=".bil$",full.names=TRUE)
prcp_files_rg <- regexpr("_([^_]*_[^_]*)$",prcp_files)
prcp_files_date <- regmatches(prcp_files,prcp_files_rg)
prcp_files_order <- order(prcp_files_date)
prcp_files <- prcp_files[prcp_files_order]
prcp_daily_avgs <- c()
for (i in 1:length(prcp_files)){
  flush.console()
  print(paste("Now processing ",prcp_files[i]))
  rast <- raster(prcp_files[i])
  avg <- cellStats(rast,'mean')
  prcp_daily_avgs <- c(prcp_daily_avgs,avg)
}

setwd("../PRISM_temp_cropped/")
temp_files <- list.files(".",pattern=".bil$",full.names=TRUE)
temp_files_rg <- regexpr("_([^_]*_[^_]*)$",temp_files)
temp_files_date <- regmatches(temp_files,temp_files_rg)
temp_files_order <- order(temp_files_date)
temp_files <- temp_files[temp_files_order]
temp_daily_avgs <- c()
for (i in 1:length(temp_files)){
  flush.console()
  print(paste("Now processing ",temp_files[i]))
  rast <- raster(temp_files[i])
  avg <- cellStats(rast,'mean')
  temp_daily_avgs <- c(temp_daily_avgs,avg)
}

setwd("../PRISM_tmin_cropped/")
tmin_files <- list.files(".",pattern=".bil$",full.names=TRUE)
tmin_files_rg <- regexpr("_([^_]*_[^_]*)$",tmin_files)
tmin_files_date <- regmatches(tmin_files,tmin_files_rg)
tmin_files_order <- order(tmin_files_date)
tmin_files <- tmin_files[tmin_files_order]
tmin_daily_avgs <- rep(NA,length(tmin_files))
for (i in 1:length(tmin_files)){
  flush.console()
  print(paste("Now processing ",tmin_files[i]))
  rast <- raster(tmin_files[i])
  tmin_daily_avgs[i] <- cellStats(rast,'mean')
}

setwd("../PRISM_tmax_cropped/")
tmax_files <- list.files(".",pattern=".bil$",full.names=TRUE)
tmax_files_rg <- regexpr("_([^_]*_[^_]*)$",tmax_files)
tmax_files_date <- regmatches(tmax_files,tmax_files_rg)
tmax_files_order <- order(tmax_files_date)
tmax_files <- tmax_files[tmax_files_order]
tmax_daily_avgs <- rep(NA,length(tmax_files))
for (i in 1:length(tmax_files)){
  flush.console()
  print(paste("Now processing ",tmax_files[i]))
  rast <- raster(tmax_files[i])
  tmax_daily_avgs[i] <- cellStats(rast,'mean')
}

## Writes to disk.
days <- seq(as.Date("2004-01-01"),as.Date("2015-6-30"),by="day")
PRISM_daily <- data.frame(DATE=days,TEMP=temp_daily_avgs,PRCP=prcp_daily_avgs,
                          TMIN=tmin_daily_avgs,TMAX=tmax_daily_avgs)
write.csv(PRISM_daily,"../PRISM_daily_2004_2015.csv")

## Plots the temp data.
PRISM_daily$TEMP <- xts(PRISM_daily$TEMP,order.by=PRISM_daily$DATE)
plot(PRISM_daily$TEMP,xlim=as.POSIXct(c("2004-1-1","2015-6-1")))

