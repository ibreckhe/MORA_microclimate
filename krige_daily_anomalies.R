##Script to estimate daily temperature anomalies using Kriging.

####Sets up workspace and loads data####
library(raster)
library(rgeos)
library(gstat)
library(ggplot2)

setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
meta <- read.csv("/Users/ian/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/airtemp_sensor_metadata_2016.csv")
meta$alt_code <- gsub(pattern="-",replacement=".",x=meta$combined_1,fixed=TRUE)

##Adds a tiny bit of noise to coordinates to prevent singularities.
meta$X_UTM <- meta$X_UTM + runif(length(meta$X_UTM),min=-0.5,max=0.5)
meta$Y_UTM <- meta$Y_UTM + runif(length(meta$Y_UTM),min=-0.5,max=0.5)
meta$utmx_s <- meta$X_UTM / 1000
meta$utmy_s <- meta$Y_UTM / 1000


source("~/code/MORA_microclimate/air_temp_functions.R")

setwd("../cleaned/")
tavg_all <- read.csv("tavg_daily_all_stations_2009_2015.csv")
tavg_all$DATE <- as.Date(tavg_all$DATE)
tmax_all <- read.csv("tmax_daily_all_stations_2009_2015.csv")
tmax_all$DATE <- as.Date(tmax_all$DATE)
tmin_all <- read.csv("tmin_daily_all_stations_2009_2015.csv")
tmin_all$DATE <- as.Date(tmin_all$DATE)

##Drops holdout testing stations.
stations_test <- c("AB08.STR.A1","AE10.A1","AM16.A1","AV14.A1","AX15.A2",
                   "HIGH.A1","LP19","OHAN.STR.A1","PA.1570m.A1","PP17.A1","SPRY.STR.A2",
                   "SUNR.STR.A1","TB13.STR.A2","TO04.A1","TO04.A2",
                   "UpperPalisades","VATR.STR.A2","VATR.UPL.A2","WHRI.UPL.A2",
                   "WON8.A1","WON9.A1")
tavg_test <- tavg_all[,c("DATE",stations_test)]
tmin_test <- tmin_all[,c("DATE",stations_test)]
tmax_test <- tmax_all[,c("DATE",stations_test)]

tavg <- tavg_all[,which(!(colnames(tavg_all) %in% stations_test))]
tmin <- tmin_all[,which(!(colnames(tmin_all) %in% stations_test))]
tmax <- tmax_all[,which(!(colnames(tmax_all) %in% stations_test))]


# tavg <- read.csv("tavg_anomalies_2009_2015.csv")
# tavg$DATE <- as.Date(tavg$date)
# tmax <- read.csv("tmax_anomalies_2009_2015.csv")
# tmax$DATE <- as.Date(tmax$date)
# tmin <- read.csv("tmin_anomalies_2009_2015.csv")
# tmin$DATE <- as.Date(tmin$date)

##Drops sites without metadata.
meta_sites <- meta$alt_code
drop_cols <- which(!(colnames(tavg) %in% meta_sites))[-1]
tavg <- tavg[,-drop_cols]
tmin <- tmin[,-drop_cols]
tmax <- tmax[,-drop_cols]

##Reads in rasters of predictors
setwd("/Volumes/ib_working/GIS/")
elev3m <- raster("mcpred_elev_3m.tif")
elev <- aggregate(elev3m,fact=30)
coldair3m <-raster("mcpred_reg_coldair_3m.tif")
coldair <- aggregate(coldair3m,fact=30) / (elev/1000)
utmx3m <- raster("mcpred_utmx_3m.tif")
utmx <- aggregate(utmx3m,fact=30)
utmy3m <- raster("mcpred_utmy_3m.tif")
utmy <- aggregate(utmy3m,fact=30)
preds <- stack(elev,coldair,utmx,utmy)
names(preds) <- c("elev","coldair","utmx","utmy")
preds$utmx_s <- preds$utmx / 1000
preds$utmy_s <- preds$utmy / 1000
pred_crs <- crs(elev3m)

##Loops through each day and produces kriging predictions.
brickpath <- "/Volumes/ib_working/MORAclim_estimates/"

##Excludes high elevations where there are no sensors.
elev_low <- preds$elev
elev_low[elev_low[]>4500] <- NA
preds$elev <- elev_low
setwd(brickpath)

##Removes tmax estimates from problem site PP17.
tmax$PP17.A1 <- NA
tmax$PP17.A2 <- NA

##Computes anomalies and daily prediction grids.
library(foreach)
library(doParallel)

#cl <- makeCluster(3)
#registerDoParallel(cl)
# foreach(i=2015:2015) %do% {
#   date_start <- as.Date(paste(i,"-01-01",sep=""))
#   date_end <- as.Date(paste(i,"-09-30",sep=""))
#   num_start <- as.numeric(format(date_start,format="%j"))
#   num_end <- as.numeric(format(date_end,format="%j"))
#   
#   tmax_exp_year <- paste("multiband_tmax_est_year_",i,".grd",sep="")
#   tmin_exp_year <- paste("multiband_tmin_est_year_",i,".grd",sep="")
#   tavg_exp_year <- paste("multiband_tavg_est_year_",i,".grd",sep="")
#   
#   tavg_year <- tavg[tavg$DATE %in% seq(date_start,date_end,by="day"),]
#   tmin_year <- tmin[tmin$DATE %in% seq(date_start,date_end,by="day"),]
#   tmax_year <- tmax[tmax$DATE %in% seq(date_start,date_end,by="day"),]
#   
#   out_string_tmax <- paste("../MORAclim_krige/krige_tmax_year_",i,sep="")
#   krige_daily_anoms(exp_brick=tmax_exp_year,exp_brick_path=brickpath,
#                     meas_df=tmax_year,
#                     meta=meta,preds=preds,
#                     out_prefix=out_string_tmax,
#                     overwrite=TRUE) 
#   out_string_tmin <- paste("../MORAclim_krige/krige_tmin_year_",i,sep="")
#   krige_daily_anoms(exp_brick=tmin_exp_year,exp_brick_path=brickpath,
#                     meas_df=tmin_year,
#                     meta=meta,preds=preds,
#                     out_prefix=out_string_tmin,
#                     overwrite=TRUE)
#   out_string_tavg <- paste("../MORAclim_krige/krige_tavg_year_",i,sep="")
#   krige_daily_anoms(exp_brick=tavg_exp_year,exp_brick_path=brickpath,
#                     meas_df=tavg_year,
#                     meta=meta,preds=preds,
#                     out_prefix=out_string_tavg,
#                     overwrite=TRUE)
#   
# }
#stopCluster(cl)

##Tries a different approach for days with large anomalies.
probs <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/Krige_test_err_days.csv")
probs$tmin_probs <- probs$avg_err_tmin > 1.5
probs$tmax_probs <- probs$avg_err_tmax > 1.5
probs$tavg_probs <- probs$avg_err_tavg > 1.5
probs$date <- as.Date(probs$date)
probs$year <- as.numeric(strftime(probs$date,format="%Y"))
probs$doy <- as.numeric(strftime(probs$date,format="%j"))

cl <- makeCluster(3)
registerDoParallel(cl)
foreach(i=2010:2014) %do% {
  date_start <- as.Date(paste(i,"-01-01",sep=""))
  date_end <- as.Date(paste(i,"-12-31",sep=""))
  num_start <- as.numeric(format(date_start,format="%j"))
  num_end <- as.numeric(format(date_end,format="%j"))

  tmax_exp_year <- paste("multiband_tmax_est_year_",i,".grd",sep="")
  tmin_exp_year <- paste("multiband_tmin_est_year_",i,".grd",sep="")
  tavg_exp_year <- paste("multiband_tavg_est_year_",i,".grd",sep="")

  tavg_year <- tavg[tavg$DATE %in% seq(date_start,date_end,by="day"),]
  tmin_year <- tmin[tmin$DATE %in% seq(date_start,date_end,by="day"),]
  tmax_year <- tmax[tmax$DATE %in% seq(date_start,date_end,by="day"),]

  probs_year <- probs[probs$date %in% seq(date_start,date_end,by="day"),]
  probs_tavg_year <- probs_year$tavg_probs
  probs_tmax_year <- probs_year$tmax_probs
  probs_tmin_year <- probs_year$tmin_probs
  
  
  out_string_tmax <- paste("../MORAclim_krige/krige_tmax_year_",i,sep="")
  krige_problem_anoms(exp_brick=tmax_exp_year,exp_brick_path=brickpath,
                    meas_df=tmax_year,probs=probs_tmax_year,
                    meta=meta,preds=preds,
                    out_prefix=out_string_tmax,
                    overwrite=TRUE)
  out_string_tmin <- paste("../MORAclim_krige/krige_tmin_year_",i,sep="")
  krige_problem_anoms(exp_brick=tmin_exp_year,exp_brick_path=brickpath,
                    meas_df=tmin_year,probs=probs_tmin_year,
                    meta=meta,preds=preds,
                    out_prefix=out_string_tmin,
                    overwrite=TRUE)
  out_string_tavg <- paste("../MORAclim_krige/krige_tavg_year_",i,sep="")
  krige_problem_anoms(exp_brick=tavg_exp_year,exp_brick_path=brickpath,
                    meas_df=tavg_year,probs=probs_tavg_year,
                    meta=meta,preds=preds,
                    out_prefix=out_string_tavg,
                    overwrite=TRUE)

}
stopCluster(cl)


##Loops through each day and produces Bayesian predictions.
setwd("/Volumes/ib_working/MORAclim_estimates/")

##Aggregates predictor grids to 180m
preds_180m <- aggregate(preds,fact=3,na.rm=TRUE)

##Computes anomalies and daily prediction grids.
cl <- makeCluster(2)
registerDoParallel(cl)
foreach(i=2015:2015) %do% {
  date_start <- as.Date(paste(i,"-01-01",sep=""))
  date_end <- as.Date(paste(i,"-09-30",sep=""))
  
  tmax_exp_year <- paste("multiband_tmax_est_year_",i,".grd",sep="")
  tmin_exp_year <- paste("multiband_tmin_est_year_",i,".grd",sep="")
  tavg_exp_year <- paste("multiband_tavg_est_year_",i,".grd",sep="")
  
  tavg_year <- tavg[tavg$DATE %in% seq(date_start,date_end,by="day"),]
  tmin_year <- tmin[tmin$DATE %in% seq(date_start,date_end,by="day"),]
  tmax_year <- tmax[tmax$DATE %in% seq(date_start,date_end,by="day"),]
  
  out_string_tmax <- paste("../MORAclim_bayes/bayes_tmax_year_",i,sep="")
  bayes_daily_anoms(exp_brick=tmax_exp_year,
                    exp_brick_path=brickpath,
                    meas_df=tmax_year,pp=FALSE,
                    meta=meta,preds=preds_180m,nsamples=1000,
                    out_prefix=out_string_tmax,
                    overwrite=TRUE)
  out_string_tmin <- paste("../MORAclim_bayes/bayes_tmin_year_",i,sep="")
  bayes_daily_anoms(exp_brick=tmin_exp_year,
                    exp_brick_path=brickpath,
                    meas_df=tmin_year,pp=FALSE,
                    meta=meta,preds=preds_180m,nsamples=1000,
                    out_prefix=out_string_tmin,
                    overwrite=TRUE)
  out_string_tavg <- paste("../MORAclim_bayes/bayes_tavg_year_",i,sep="")
  bayes_daily_anoms(exp_brick=tavg_exp_year,
                    exp_brick_path=brickpath,
                    meas_df=tavg_year,pp=FALSE,
                    meta=meta,preds=preds_180m,nsamples=1000,
                    out_prefix=out_string_tavg,
                    overwrite=TRUE)
}
stopCluster(cl)

####Puts all files in yearly raster bricks.####



