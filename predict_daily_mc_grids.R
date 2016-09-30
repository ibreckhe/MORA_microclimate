##Script to develop daily grids of expected microclimate.
##Author: Ian Breckheimer
##Date: 9 February 2015

##Sets up workspace
library(raster)
library(xts)

####Brings in data.####
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
source("~/code/MORA_microclimate/air_temp_functions.R")

##PRISM series
prism <- read.csv("PRISM_daily_1_2004_9_2015.csv")
prism$DATE <- as.Date(prism$DATE)

##PRISM Seasonal Averages.
prism_interp <- read.csv("prism_seas_avg_interp_2009_2015.csv")
prism_interp$date <- as.Date(prism_interp$date)

##day of year estimates.
prism_year_avg <- prism_interp[prism_interp$date %in% seq(as.Date("2010-01-01"),as.Date("2010-12-31"),by="day"),]

##Daily temperature anomalies.
tavg_anoms <- read.csv("../cleaned/tavg_anomalies_2009_2015.csv")
tmin_anoms <- read.csv("../cleaned/tmin_anomalies_2009_2015.csv")
tmax_anoms <- read.csv("../cleaned/tmax_anomalies_2009_2015.csv")

##Brings in test raster data.
#test_disp <- brick("~/GIS/test_disp.grd")
#test_sens <- brick("~/GIS/test_sens.grd")

# ##Subsets estimates for each climate variable.
# test_disp_tavg <- test_disp[[c(5,11,17,23)]]
# test_disp_tmin <- test_disp[[c(3,9,15,21)]]
# test_disp_tmax <- test_disp[[c(1,7,13,19)]]
# 
# test_sens_tavg <- test_sens[[c(5,11,17,23)]]
# test_sens_tmin <- test_sens[[c(3,9,15,21)]]
# test_sens_tmax <- test_sens[[c(1,7,13,19)]]

##Brings in full extent rasters.

names <- c("DJF","MAM","JJA","SON")

##Disparity
setwd("/Volumes/ib_working/GIS/")
disp_tavg_files <- list.files(".",pattern="disp_tavg_gam.grd$")
disp_tavg <- stack(disp_tavg_files)
disp_tavg <- disp_tavg[[c(1,5,3,7)]]
names(disp_tavg) <- paste("disp_tavg_",names,sep="")

disp_tmax_files <- list.files(".",pattern="disp_tmax_gam.grd$")
disp_tmax <- stack(disp_tmax_files)
disp_tmax <- disp_tmax[[c(1,5,3,7)]]
names(disp_tmax) <- paste("disp_tmax_",names,sep="")

disp_tmin_files <- list.files(".",pattern="disp_tmin_gam.grd$")
disp_tmin <- stack(disp_tmin_files)
disp_tmin <- disp_tmin[[c(1,5,3,7)]]
names(disp_tmin) <- paste("disp_tmin_",names,sep="")

##Sensitivity.
sens_tavg_files <- list.files(".",pattern="sens_tavg_gam.grd$")
sens_tavg <- stack(sens_tavg_files)
sens_tavg <- sens_tavg[[c(1,5,3,7)]]
names(sens_tavg) <- paste("sens_tavg_",names,sep="")

sens_tmax_files <- list.files(".",pattern="sens_tmax_gam.grd$")
sens_tmax <- stack(sens_tmax_files)
sens_tmax <- sens_tmax[[c(1,5,3,7)]]
names(sens_tmax) <- paste("sens_tmax_",names,sep="")

sens_tmin_files <- list.files(".",pattern="sens_tmin_gam.grd$")
sens_tmin <- stack(sens_tmin_files)
sens_tmin <- sens_tmin[[c(1,5,3,7)]]
names(sens_tmin) <- paste("sens_tmin_",names,sep="")

##Aggregates them to 90m resolution.
disp_tmax_90m <- aggregate(disp_tmax,fact=30,expand=TRUE,overwrite=TRUE,fun=mean,
                           na.rm=TRUE,filename="~/GIS/clim_seas_disp_tmax_gam_90m.grd")
disp_tmin_90m <- aggregate(disp_tmin,fact=30,expand=TRUE,overwrite=TRUE,fun=mean,
                           na.rm=TRUE,filename="~/GIS/clim_seas_disp_tmin_gam_90m.grd")
disp_tavg_90m <- aggregate(disp_tavg,fact=30,expand=TRUE,overwrite=TRUE,fun=mean,
                           na.rm=TRUE,filename="~/GIS/clim_seas_disp_tavg_gam_90m.grd")
sens_tmax_90m <- aggregate(sens_tmax,fact=30,expand=TRUE,overwrite=TRUE,fun=mean,
                           na.rm=TRUE,filename="~/GIS/clim_seas_sens_tmax_gam_90m.grd")
sens_tmin_90m <- aggregate(sens_tmin,fact=30,expand=TRUE,overwrite=TRUE,fun=mean,
                           na.rm=TRUE,filename="~/GIS/clim_seas_sens_tmin_gam_90m.grd")
sens_tavg_90m <- aggregate(sens_tavg,fact=30,expand=TRUE,overwrite=TRUE,fun=mean,
                           na.rm=TRUE,filename="~/GIS/clim_seas_sens_tavg_gam_90m.grd")

####Linear interpolation of model coefficients for each day.####

##Sequence of days starting on October 15th and ending January 15th of next year

##Preps stack for interpolation.
disp_tmax_90m_int <- disp_tmax_90m[[c("disp_tmax_SON",
                                  "disp_tmax_DJF",
                                  "disp_tmax_MAM",
                                  "disp_tmax_JJA",
                                  "disp_tmax_SON",
                                  "disp_tmax_DJF")]]

disp_tmin_90m_int <- disp_tmin_90m[[c("disp_tmin_SON",
                                      "disp_tmin_DJF",
                                      "disp_tmin_MAM",
                                      "disp_tmin_JJA",
                                      "disp_tmin_SON",
                                      "disp_tmin_DJF")]]

disp_tavg_90m_int <- disp_tavg_90m[[c("disp_tavg_SON",
                                      "disp_tavg_DJF",
                                      "disp_tavg_MAM",
                                      "disp_tavg_JJA",
                                      "disp_tavg_SON",
                                      "disp_tavg_DJF")]]

sens_tmax_90m_int <- sens_tmax_90m[[c("sens_tmax_SON",
                                      "sens_tmax_DJF",
                                      "sens_tmax_MAM",
                                      "sens_tmax_JJA",
                                      "sens_tmax_SON",
                                      "sens_tmax_DJF")]]

sens_tmin_90m_int <- sens_tmin_90m[[c("sens_tmin_SON",
                                      "sens_tmin_DJF",
                                      "sens_tmin_MAM",
                                      "sens_tmin_JJA",
                                      "sens_tmin_SON",
                                      "sens_tmin_DJF")]]

sens_tavg_90m_int <- sens_tavg_90m[[c("sens_tavg_SON",
                                      "sens_tavg_DJF",
                                      "sens_tavg_MAM",
                                      "sens_tavg_JJA",
                                      "sens_tavg_SON",
                                      "sens_tavg_DJF")]]

day_dates <- seq(as.Date("2013-10-15"),as.Date("2015-01-15"),by="day")
day_seq <- 1:length(day_dates)
seas_days <- c(1,93,183,274,366,458)

##Interpolates the rasters

disp_tmax_interp <- interpolateTemporal(disp_tmax_90m_int,xin=seas_days,xout=day_seq,
                                   prefix="disp_tmax_interp",writechange=FALSE,overwrite=TRUE,
                                   outdir="/Volumes/ib_working/GIS/microclim_interp/")
disp_tmin_interp <- interpolateTemporal(disp_tmin_90m_int,xin=seas_days,xout=day_seq,
                                        prefix="disp_tmin_interp",writechange=FALSE,overwrite=TRUE,
                                        outdir="/Volumes/ib_working/GIS/microclim_interp/")
disp_tavg_interp <- interpolateTemporal(disp_tavg_90m_int,xin=seas_days,xout=day_seq,
                                        prefix="disp_tavg_interp",writechange=FALSE,overwrite=TRUE,
                                        outdir="/Volumes/ib_working/GIS/microclim_interp/")

sens_tmax_interp <- interpolateTemporal(sens_tmax_90m_int,xin=seas_days,xout=day_seq,
                                        prefix="sens_tmax_interp",writechange=FALSE,overwrite=TRUE,
                                        outdir="/Volumes/ib_working/GIS/microclim_interp/")
sens_tmin_interp <- interpolateTemporal(sens_tmin_90m_int,xin=seas_days,xout=day_seq,
                                        prefix="sens_tmin_interp",writechange=FALSE,overwrite=TRUE,
                                        outdir="/Volumes/ib_working/GIS/microclim_interp/")
sens_tavg_interp <- interpolateTemporal(sens_tavg_90m_int,xin=seas_days,xout=day_seq,
                                        prefix="sens_tavg_interp",writechange=FALSE,overwrite=TRUE,
                                        outdir="/Volumes/ib_working/GIS/microclim_interp/")

##Brings interpolated grids back in as a stack.
setwd("/Volumes/ib_working/GIS/microclim_interp/")
disp_tmax_files <- list.files(".",pattern="disp_tmax_interp")
disp_tmax_interp <- stack(disp_tmax_files,quick=TRUE)
disp_tmin_files <- list.files(".",pattern="disp_tmin_interp")
disp_tmin_interp <- stack(disp_tmin_files,quick=TRUE)
disp_tavg_files <- list.files(".",pattern="disp_tavg_interp")
disp_tavg_interp <- stack(disp_tavg_files,quick=TRUE)

sens_tmax_files <- list.files(".",pattern="sens_tmax_interp")
sens_tmax_interp <- stack(sens_tmax_files,quick=TRUE)
sens_tmin_files <- list.files(".",pattern="sens_tmin_interp")
sens_tmin_interp <- stack(sens_tmin_files,quick=TRUE)
sens_tavg_files <- list.files(".",pattern="sens_tavg_interp")
sens_tavg_interp <- stack(sens_tavg_files,quick=TRUE)

##Subsets for just January 1 to December 31.
daynames <- paste("day_",1:365,sep="")
dayindex <- 79:443

disp_tmax_interp <- disp_tmax_interp[[dayindex]]
names(disp_tmax_interp) <- daynames
disp_tmin_interp <- disp_tmin_interp[[dayindex]]
names(disp_tmin_interp) <- daynames
disp_tavg_interp <- disp_tavg_interp[[dayindex]]
names(disp_tavg_interp) <- daynames

sens_tmax_interp <- sens_tmax_interp[[dayindex]]
names(sens_tmax_interp) <- daynames
sens_tmin_interp <- sens_tmin_interp[[dayindex]]
names(sens_tmin_interp) <- daynames
sens_tavg_interp <- sens_tavg_interp[[dayindex]]
names(sens_tavg_interp) <- daynames

##Writes raster bricks of interpolated coefficients to disk.
disp_tmax_interp_brick <- brick(disp_tmax_interp,overwrite=TRUE,
                                filename="../interp_disp_tmax_90m.grd")
disp_tmin_interp_brick <- brick(disp_tmin_interp,overwrite=TRUE,
                                filename="../interp_disp_tmin_90m.grd")
disp_tavg_interp_brick <- brick(disp_tavg_interp,overwrite=TRUE,
                                filename="../interp_disp_tavg_90m.grd")

sens_tmax_interp_brick <- brick(sens_tmax_interp,overwrite=TRUE,
                                filename="../interp_sens_tmax_90m.grd")
sens_tmin_interp_brick <- brick(sens_tmin_interp,overwrite=TRUE,
                                filename="../interp_sens_tmin_90m.grd")
sens_tavg_interp_brick <- brick(sens_tavg_interp,overwrite=TRUE,
                                filename="../interp_sens_tavg_90m.grd")

##Reads them back in
setwd("/Volumes/ib_working/GIS/microclim_interp/")
disp_tmax_interp_brick <- brick("../interp_disp_tmax_90m.grd")
disp_tmin_interp_brick <- brick("../interp_disp_tmin_90m.grd")
disp_tavg_interp_brick <- brick("../interp_disp_tavg_90m.grd")

sens_tmax_interp_brick <- brick("../interp_sens_tmax_90m.grd")
sens_tmin_interp_brick <- brick("../interp_sens_tmin_90m.grd")
sens_tavg_interp_brick <- brick("../interp_sens_tavg_90m.grd")

####Creates grids of expected climate variables for each day from 2009 to 2015.####

##Creates raster stacks with constants for each day.
prism_tavg_leap <- disp_tavg_interp_brick
prism_tavg_leap$day_366 <- NA
prism_tmin_leap <- disp_tavg_interp_brick
prism_tmin_leap$day_366 <- NA
prism_tmax_leap <- disp_tavg_interp_brick
prism_tmax_leap$day_366 <- NA
tmax_leap <- c(prism_year_avg$tmax_seas,last(prism_year_avg$tmax_seas))
tmin_leap <- c(prism_year_avg$tmin_seas,last(prism_year_avg$tmin_seas))
tavg_leap <- c(prism_year_avg$tavg_seas,last(prism_year_avg$tavg_seas))

disp_tmax_interp_leap <- disp_tmax_interp_brick
disp_tmax_interp_leap$day_366 <- disp_tmax_interp_brick[[365]]
disp_tmin_interp_leap <- disp_tmin_interp_brick
disp_tmin_interp_leap$day_366 <- disp_tmin_interp_brick[[365]]
disp_tavg_interp_leap <- disp_tavg_interp_brick
disp_tavg_interp_leap$day_366 <- disp_tavg_interp_brick[[365]]

sens_tmax_interp_leap <- sens_tmax_interp_brick
sens_tmax_interp_leap$day_366 <- sens_tmax_interp_brick[[365]]
sens_tmin_interp_leap <- sens_tmin_interp_brick
sens_tmin_interp_leap$day_366 <- sens_tmin_interp_brick[[365]]
sens_tavg_interp_leap <- sens_tavg_interp_brick
sens_tavg_interp_leap$day_366 <- sens_tavg_interp_brick[[365]]


####Creates gridded 90m expectations for each day from 2009 to 2015.####
reg_fun <- function(a,b,c,d) {((a - b) * c + d + b) * 10000}
setwd("/Volumes/ib_working/MORAclim_estimates/")
for (j in c(2015:2015)){
  print(paste("Processing year",j))
  
  start_date <- as.Date(paste(j,"-01-01",sep=""))
  end_date <- as.Date(paste(j,"-12-31",sep=""))
  reg_date <- prism[prism$DATE %in% seq(start_date,end_date,by="day"),]
  
  ##Leap year detector
  if(j/4 == round(j/4)){
    for (i in 1:366){
      print(paste("Leap year, processing day ",i," of ",366,sep=""))
      
      prism_rast_tmax <-  prism_tmax_leap[[i]]
      prism_rast_tmax[] <- reg_date$TMAX[i]
      prism_avg_tmax <- prism_tmax_leap[[i]]
      prism_avg_tmax[] <- tmax_leap[i]
      filename <- paste("tmax_est_year_",j,"_day_",sprintf("%04d", i),".tif",sep="")
      out_est_tmax <- overlay(prism_rast_tmax,prism_avg_tmax,sens_tmax_interp_leap[[i]],
                              disp_tmax_interp_leap[[i]],fun=reg_fun,datatype="INT4S",
                              filename=filename,unstack=FALSE,overwrite=TRUE)
      
      prism_rast_tmin <-  prism_tmin_leap[[i]]
      prism_rast_tmin[] <- reg_date$TMIN[i]
      prism_avg_tmin <- prism_tmin_leap[[i]]
      prism_avg_tmin[] <- tmin_leap[i]
      filename <- paste("tmin_est_year_",j,"_day_",sprintf("%04d", i),".tif",sep="")
      out_est_tmin <- overlay(prism_rast_tmin,prism_avg_tmin,sens_tmin_interp_leap[[i]],
                              disp_tmin_interp_leap[[i]],fun=reg_fun,datatype="INT4S",
                              filename=filename,unstack=FALSE,overwrite=TRUE)
      
      prism_rast_tavg <-  prism_tavg_leap[[i]]
      prism_rast_tavg[] <- reg_date$TEMP[i]
      prism_avg_tavg <- prism_tavg_leap[[i]]
      prism_avg_tavg[] <- tavg_leap[i]
      filename <- paste("tavg_est_year_",j,"_day_",sprintf("%04d", i),".tif",sep="")
      out_est_tavg <- overlay(prism_rast_tavg,prism_avg_tavg,sens_tavg_interp_leap[[i]],
                              disp_tavg_interp_leap[[i]],fun=reg_fun,datatype="INT4S",
                              filename=filename,unstack=FALSE,overwrite=TRUE)
    }
  }else{
    for(i in 1:365){
      print(paste("Processing day ",i," of ",nlayers(disp_tmax_interp_brick),sep=""))
      
      prism_rast_tmax <-  disp_tmax_interp_brick[[i]]
      prism_rast_tmax[] <- reg_date$TMAX[i]
      prism_avg_tmax <- disp_tmax_interp_brick[[i]]
      prism_avg_tmax[] <- prism_year_avg$tmax_seas[i]
      filename <- paste("tmax_est_year_",j,"_day_",sprintf("%04d", i),".tif",sep="")
      out_est_tmax <- overlay(prism_rast_tmax,prism_avg_tmax,sens_tmax_interp_brick[[i]],
                              disp_tmax_interp_brick[[i]],fun=reg_fun,datatype="INT4S",
                              filename=filename,unstack=FALSE,overwrite=TRUE)
      
      prism_rast_tmin <-  disp_tmin_interp_brick[[i]]
      prism_rast_tmin[] <- reg_date$TMIN[i]
      prism_avg_tmin <- disp_tmin_interp_brick[[i]]
      prism_avg_tmin[] <- prism_year_avg$tmin_seas[i]
      filename <- paste("tmin_est_year_",j,"_day_",sprintf("%04d", i),".tif",sep="")
      out_est_tmin <- overlay(prism_rast_tmin,prism_avg_tmin,sens_tmin_interp_brick[[i]],
                              disp_tmin_interp_brick[[i]],fun=reg_fun,datatype="INT4S",
                              filename=filename,unstack=FALSE,overwrite=TRUE)
      
      prism_rast_tavg <-  disp_tavg_interp_brick[[i]]
      prism_rast_tavg[] <- reg_date$TEMP[i]
      prism_avg_tavg <- disp_tavg_interp_brick[[i]]
      prism_avg_tavg[] <- prism_year_avg$tavg_seas[i]
      filename <- paste("tavg_est_year_",j,"_day_",sprintf("%04d", i),".tif",sep="")
      out_est_tavg <- overlay(prism_rast_tavg,prism_avg_tavg,sens_tavg_interp_brick[[i]],
                              disp_tavg_interp_brick[[i]],fun=reg_fun,datatype="INT4S",
                              filename=filename,unstack=FALSE,overwrite=TRUE)
    }
  }
}

##Puts data for each year in a raster brick.
for (k in c(2009:2015)){
  year_pat_tavg <- paste("tavg_est_year_",k,sep="")
  year_list_tavg <- list.files(".",pattern=year_pat_tavg)
  year_stack_tavg <- stack(year_list_tavg,quick=TRUE)
  year_brick_tavg <- brick(year_stack_tavg,overwrite=TRUE,
                      filename=paste("multiband_",year_pat_tavg,".grd",sep=""))
  
  year_pat_tmax <- paste("tmax_est_year_",k,sep="")
  year_list_tmax <- list.files(".",pattern=year_pat_tmax)
  year_stack_tmax <- stack(year_list_tmax,quick=TRUE)
  year_brick_tmax <- brick(year_stack_tmax,overwrite=TRUE,
                           filename=paste("multiband_",year_pat_tmax,".grd",sep=""))
  year_pat_tmin <- paste("tmin_est_year_",k,sep="")
  year_list_tmin <- list.files(".",pattern=year_pat_tmin)
  year_stack_tmin <- stack(year_list_tmin,quick=TRUE)
  year_brick_tmin <- brick(year_stack_tmin,overwrite=TRUE,
                           filename=paste("multiband_",year_pat_tmin,".grd",sep=""))
}

##Extracts data to check predictions against observed data.

tavg_estimates <- brick("multiband_tavg_est_year_2014.grd")


setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/cleaned_dailyair/")
tavg <- read.csv("alldat_daily_tavg_2006_2015.csv")
tavg$DATE <- as.Date(tavg$DATE)
tmax <- read.csv("alldat_daily_tmax_2006_2015.csv")
tmax$DATE <- as.Date(tmax$DATE)
tmin <- read.csv("alldat_daily_tmin_2006_2015.csv")
tmin$DATE <- as.Date(tmin$DATE)
meta <- read.csv("/Users/ian/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/airtemp_site_metadata_2015.csv")
meta$alt_code <- gsub(pattern="-",replacement=".",x=meta$location,fixed=TRUE)

test_coords <- c(meta$utmx[meta$alt_code=="TO04.A1"],meta$utmy[meta$alt_code=="TO04.A1"])
test_reg <- as.numeric(extract(tavg_estimates,cbind(test_coords[1],test_coords[2]))) / 10000
tavg_2014 <- tavg[tavg$DATE %in% seq(as.Date("2014-01-01"),as.Date("2014-12-31"),by="day"),]
test_real <- tavg_2014$TO04.A1

plot(test_reg[-1],test_real[1:(length(test_real)-1)])
abline(c(0,1),lty=2)
