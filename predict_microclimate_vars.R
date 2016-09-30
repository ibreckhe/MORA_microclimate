##Script to predict rasters of microclimate variables.
##Author:Ian Breckheimer
##7 February 2015.

library(mgcv)
library(ggplot2)
library(gridExtra)
library(raster)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw")

##Reads from disk.
meta_indices_tmax <- read.csv("../cleaned/meta_indices_dry_tmax_2009_2015.csv")
meta_indices_tmin <- read.csv("../cleaned/meta_indices_wet_tmin_2009_2015.csv")
meta_indices_tavg <- read.csv("../cleaned/meta_indices_tavg_2009_2015.csv")

##Reads seasonal estimates from disk.
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

##Joins with covariates.
meta_indices_DJF_tavg <- merge(meta_indices_tavg[,1:42],indices_DJF_tavg,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_DJF_tmin <- merge(meta_indices_tmin[,1:42],indices_DJF_tmin,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_DJF_tmax <- merge(meta_indices_tmax[,1:42],indices_DJF_tmax,by.x="alt_code",by.y="site",keep.x=TRUE)

meta_indices_MAM_tavg <- merge(meta_indices_tavg[,1:42],indices_MAM_tavg,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_MAM_tmin <- merge(meta_indices_tmin[,1:42],indices_MAM_tmin,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_MAM_tmax <- merge(meta_indices_tmax[,1:42],indices_MAM_tmax,by.x="alt_code",by.y="site",keep.x=TRUE)

meta_indices_JJA_tavg <- merge(meta_indices_tavg[,1:42],indices_JJA_tavg,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_JJA_tmin <- merge(meta_indices_tmin[,1:42],indices_JJA_tmin,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_JJA_tmax <- merge(meta_indices_tmax[,1:42],indices_JJA_tmax,by.x="alt_code",by.y="site",keep.x=TRUE)

meta_indices_SON_tavg <- merge(meta_indices_tavg[,1:42],indices_SON_tavg,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_SON_tmin <- merge(meta_indices_tmin[,1:42],indices_SON_tmin,by.x="alt_code",by.y="site",keep.x=TRUE)
meta_indices_SON_tmax <- merge(meta_indices_tmax[,1:42],indices_SON_tmax,by.x="alt_code",by.y="site",keep.x=TRUE)

##Drops sites with large uncertainties.
filter_thresh <- 0.8

meta_indices_tavg <- meta_indices_tavg[(meta_indices_tavg$disp_upr95 - meta_indices_tavg$disp_lwr95) < filter_thresh,]
meta_indices_tmin <- meta_indices_tmin[(meta_indices_tmin$disp_upr95 - meta_indices_tmin$disp_lwr95) < filter_thresh,]
meta_indices_tmax <- meta_indices_tmax[(meta_indices_tmax$disp_upr95 - meta_indices_tmax$disp_lwr95) < filter_thresh,]

meta_indices_DJF_tavg <- meta_indices_DJF_tavg[(meta_indices_DJF_tavg$disp_upr95 - meta_indices_DJF_tavg$disp_lwr95) < filter_thresh,]
meta_indices_DJF_tmin <- meta_indices_DJF_tmin[(meta_indices_DJF_tmin$disp_upr95 - meta_indices_DJF_tmin$disp_lwr95) < filter_thresh,]
meta_indices_DJF_tmax <- meta_indices_DJF_tmax[(meta_indices_DJF_tmax$disp_upr95 - meta_indices_DJF_tmax$disp_lwr95) < filter_thresh,]

meta_indices_MAM_tavg <- meta_indices_MAM_tavg[(meta_indices_MAM_tavg$disp_upr95 - meta_indices_MAM_tavg$disp_lwr95) < filter_thresh,]
meta_indices_MAM_tmin <- meta_indices_MAM_tmin[(meta_indices_MAM_tmin$disp_upr95 - meta_indices_MAM_tmin$disp_lwr95) < filter_thresh,]
meta_indices_MAM_tmax <- meta_indices_MAM_tmax[(meta_indices_MAM_tmax$disp_upr95 - meta_indices_MAM_tmax$disp_lwr95) < filter_thresh,]

meta_indices_JJA_tavg <- meta_indices_JJA_tavg[(meta_indices_JJA_tavg$disp_upr95 - meta_indices_JJA_tavg$disp_lwr95) < filter_thresh,]
meta_indices_JJA_tmin <- meta_indices_JJA_tmin[(meta_indices_JJA_tmin$disp_upr95 - meta_indices_JJA_tmin$disp_lwr95) < filter_thresh,]
meta_indices_JJA_tmax <- meta_indices_JJA_tmax[(meta_indices_JJA_tmax$disp_upr95 - meta_indices_JJA_tmax$disp_lwr95) < filter_thresh,]

meta_indices_SON_tavg <- meta_indices_SON_tavg[(meta_indices_SON_tavg$disp_upr95 - meta_indices_SON_tavg$disp_lwr95) < filter_thresh,]
meta_indices_SON_tmin <- meta_indices_SON_tmin[(meta_indices_SON_tmin$disp_upr95 - meta_indices_SON_tmin$disp_lwr95) < filter_thresh,]
meta_indices_SON_tmax <- meta_indices_SON_tmax[(meta_indices_SON_tmax$disp_upr95 - meta_indices_SON_tmax$disp_lwr95) < filter_thresh,]

##Drops suspect stations.
drop_names <- c("MU101-A1","PACW1-A1")

meta_indices_tmax <- meta_indices_tmax[-c(which(meta_indices_tmax$combined_n %in% drop_names)),]
meta_indices_tmin<- meta_indices_tmin[-c(which(meta_indices_tmin$combined_n %in% drop_names)),]
meta_indices_tavg <- meta_indices_tavg[-c(which(meta_indices_tavg$combined_n %in% drop_names)),]

meta_indices_DJF_tmax <- meta_indices_DJF_tmax[-c(which(meta_indices_DJF_tmax$combined_n %in% drop_names)),]
meta_indices_DJF_tmin<- meta_indices_DJF_tmin[-c(which(meta_indices_DJF_tmin$combined_n %in% drop_names)),]
meta_indices_DJF_tavg <- meta_indices_DJF_tavg[-c(which(meta_indices_DJF_tavg$combined_n %in% drop_names)),]

meta_indices_MAM_tmax <- meta_indices_MAM_tmax[-c(which(meta_indices_MAM_tmax$combined_n %in% drop_names)),]
meta_indices_MAM_tmin<- meta_indices_MAM_tmin[-c(which(meta_indices_MAM_tmin$combined_n %in% drop_names)),]
meta_indices_MAM_tavg <- meta_indices_MAM_tavg[-c(which(meta_indices_MAM_tavg$combined_n %in% drop_names)),]

meta_indices_JJA_tmax <- meta_indices_JJA_tmax[-c(which(meta_indices_JJA_tmax$combined_n %in% drop_names)),]
meta_indices_JJA_tmin<- meta_indices_JJA_tmin[-c(which(meta_indices_JJA_tmin$combined_n %in% drop_names)),]
meta_indices_JJA_tavg <- meta_indices_JJA_tavg[-c(which(meta_indices_JJA_tavg$combined_n %in% drop_names)),]

meta_indices_SON_tmax <- meta_indices_SON_tmax[-c(which(meta_indices_SON_tmax$combined_n %in% drop_names)),]
meta_indices_SON_tmin<- meta_indices_SON_tmin[-c(which(meta_indices_SON_tmin$combined_n %in% drop_names)),]
meta_indices_SON_tavg <- meta_indices_SON_tavg[-c(which(meta_indices_SON_tavg$combined_n %in% drop_names)),]

##Reserves a random sample of 20% of stations for testing.
set.seed(25)
stations <- unique(meta_indices_tavg$combined_1)
rand <- runif(length(stations),0,1)
stations_test <- as.character(stations[rand>=0.8])
# stations_test <- c("AB08-STR-A1","AE10-A1","AM16-A1","AV14-A1","AX15-A2",
#                    "HIGH-A1","LP19","OHAN-STR-A1","PA-1570m-A1","PP17-A1","SPRY-STR-A2",
#                    "SUNR-STR-A1","TB13-STR-A2","TO04-A1","TO04-A2","TWHIT",
#                    "UpperPalisades","VATR-STR-A2","VATR-UPL-A2","WHRI-UPL-A2",
#                    "WON8-A1","WON9-A1")
stations_train <- as.character(stations[rand<0.8])

meta_indices_tavg <- meta_indices_tavg[meta_indices_tavg$combined_1 %in% stations_train,]
meta_indices_tmin <- meta_indices_tmin[meta_indices_tmin$combined_1 %in% stations_train,]
meta_indices_tmax <- meta_indices_tmax[meta_indices_tmax$combined_1 %in% stations_train,]

meta_indices_DJF_tavg <- meta_indices_DJF_tavg[meta_indices_DJF_tavg$combined_1 %in% stations_train,]
meta_indices_DJF_tmin <- meta_indices_DJF_tmin[meta_indices_DJF_tmin$combined_1 %in% stations_train,]
meta_indices_DJF_tmax <- meta_indices_DJF_tmax[meta_indices_DJF_tmax$combined_1 %in% stations_train,]

meta_indices_MAM_tavg <- meta_indices_MAM_tavg[meta_indices_MAM_tavg$combined_1 %in% stations_train,]
meta_indices_MAM_tmin <- meta_indices_MAM_tmin[meta_indices_MAM_tmin$combined_1 %in% stations_train,]
meta_indices_MAM_tmax <- meta_indices_MAM_tmax[meta_indices_MAM_tmax$combined_1 %in% stations_train,]

meta_indices_JJA_tavg <- meta_indices_JJA_tavg[meta_indices_JJA_tavg$combined_1 %in% stations_train,]
meta_indices_JJA_tmin <- meta_indices_JJA_tmin[meta_indices_JJA_tmin$combined_1 %in% stations_train,]
meta_indices_JJA_tmax <- meta_indices_JJA_tmax[meta_indices_JJA_tmax$combined_1 %in% stations_train,]

meta_indices_SON_tavg <- meta_indices_SON_tavg[meta_indices_SON_tavg$combined_1 %in% stations_train,]
meta_indices_SON_tmin <- meta_indices_SON_tmin[meta_indices_SON_tmin$combined_1 %in% stations_train,]
meta_indices_SON_tmax <- meta_indices_SON_tmax[meta_indices_SON_tmax$combined_1 %in% stations_train,]

##Fits annual model for disparity and sensitivity.
disp_tavg_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=5)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_tavg)
summary(disp_tavg_gam)

sens_tavg_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=5)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                     data=meta_indices_tavg)
summary(sens_tavg_gam)

disp_wet_tmin_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=5)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                     data=meta_indices_tmin)
summary(disp_wet_tmin_gam)

sens_wet_tmin_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=5)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                     data=meta_indices_tmin)
summary(sens_wet_tmin_gam)

disp_dry_tmax_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=5)+I(srad_seas_JJA/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_tmax)
summary(disp_dry_tmax_gam)

sens_dry_tmax_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=5)+I(srad_seas_JJA/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                     data=meta_indices_tmax)
summary(sens_dry_tmax_gam)


##Fits seasonal models for disparity and sensitivity.
DJF_disp_tmax_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_DJF/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_DJF_tmax)
summary(DJF_disp_tmax_gam)
MAM_disp_tmax_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_MAM/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_MAM_tmax)
summary(MAM_disp_tmax_gam)
JJA_disp_tmax_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_JJA/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_JJA_tmax)
summary(JJA_disp_tmax_gam)
SON_disp_tmax_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_SON/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_SON_tmax)
summary(SON_disp_tmax_gam)

DJF_disp_tmin_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_DJF_tmin)
summary(DJF_disp_tmin_gam)
MAM_disp_tmin_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_MAM_tmin)
summary(MAM_disp_tmin_gam)
JJA_disp_tmin_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_JJA_tmin)
summary(JJA_disp_tmin_gam)
SON_disp_tmin_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_SON_tmin)
summary(SON_disp_tmin_gam)

DJF_disp_tavg_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_DJF/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_DJF_tavg)
summary(DJF_disp_tavg_gam)
MAM_disp_tavg_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_MAM/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_MAM_tavg)
summary(MAM_disp_tavg_gam)
JJA_disp_tavg_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_JJA/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_JJA_tavg)
summary(JJA_disp_tavg_gam)
SON_disp_tavg_gam <- gam(disp~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_SON/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_SON_tavg)
summary(SON_disp_tavg_gam)

##Fits seasonal models for sensitivity.
DJF_sens_tmax_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_DJF/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_DJF_tmax)
summary(DJF_sens_tmax_gam)
MAM_sens_tmax_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_MAM/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_MAM_tmax)
summary(MAM_sens_tmax_gam)
JJA_sens_tmax_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_JJA/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_JJA_tmax)
summary(JJA_sens_tmax_gam)
SON_sens_tmax_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_SON/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_SON_tmax)
summary(SON_sens_tmax_gam)

##Fits seasonal models for sensitivity.
DJF_sens_tmin_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_DJF_tmin)
summary(DJF_sens_tmin_gam)
MAM_sens_tmin_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_MAM_tmin)
summary(MAM_sens_tmin_gam)
JJA_sens_tmin_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_JJA_tmin)
summary(JJA_sens_tmin_gam)
SON_sens_tmin_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_SON_tmin)
summary(SON_sens_tmin_gam)

##Fits seasonal models for sensitivity.
DJF_sens_tavg_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_DJF/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_DJF_tavg)
summary(DJF_sens_tavg_gam)
MAM_sens_tavg_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_MAM/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_MAM_tavg)
summary(MAM_sens_tavg_gam)
JJA_sens_tavg_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_JJA/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_JJA_tavg)
summary(JJA_sens_tavg_gam)
SON_sens_tavg_gam <- gam(sens~s(elev,cair_reg_9m,k=5)+s(canvol_81m,k=3)+I(srad_seas_SON/1e6)+I(X_UTM/1000)+I(Y_UTM/1000)+type,
                         data=meta_indices_SON_tavg)
summary(SON_sens_tavg_gam)

##Puts all the annual disp model objects in a list.
ann_disp_gam_mods <- list(disp_tavg_gam,disp_wet_tmin_gam,disp_dry_tmax_gam)
ann_sens_gam_mods <- list(sens_tavg_gam,sens_wet_tmin_gam,sens_dry_tmax_gam)

##Puts all the seasonal disp model objects in a list.
seas_disp_gam_mods <- list(DJF_disp_tmax_gam=DJF_disp_tmax_gam,
                           DJF_disp_tmin_gam=DJF_disp_tmin_gam,
                           DJF_disp_tavg_gam=DJF_disp_tavg_gam,
                           MAM_disp_tmax_gam=MAM_disp_tmax_gam,
                           MAM_disp_tmin_gam=MAM_disp_tmin_gam,
                           MAM_disp_tavg_gam=MAM_disp_tavg_gam,
                           JJA_disp_tmax_gam=JJA_disp_tmax_gam,
                           JJA_disp_tmin_gam=JJA_disp_tmin_gam,
                           JJA_disp_tavg_gam=JJA_disp_tavg_gam,
                           SON_disp_tmax_gam=SON_disp_tmax_gam,
                           SON_disp_tmin_gam=SON_disp_tmin_gam,
                           SON_disp_tavg_gam=SON_disp_tavg_gam)

seas_sens_gam_mods <- list(DJF_sens_tmax_gam=DJF_sens_tmax_gam,
                           DJF_sens_tmin_gam=DJF_sens_tmin_gam,
                           DJF_sens_tavg_gam=DJF_sens_tavg_gam,
                           MAM_sens_tmax_gam=MAM_sens_tmax_gam,
                           MAM_sens_tmin_gam=MAM_sens_tmin_gam,
                           MAM_sens_tavg_gam=MAM_sens_tavg_gam,
                           JJA_sens_tmax_gam=JJA_sens_tmax_gam,
                           JJA_sens_tmin_gam=JJA_sens_tmin_gam,
                           JJA_sens_tavg_gam=JJA_sens_tavg_gam,
                           SON_sens_tmax_gam=SON_sens_tmax_gam,
                           SON_sens_tmin_gam=SON_sens_tmin_gam,
                           SON_sens_tavg_gam=SON_sens_tavg_gam)

clim_brick <- brick("~/GIS/moraclim_preds_2016.grd")
names(clim_brick) <- c("X_UTM","Y_UTM","canvol_81m","can_pct_81m","elev","cair_reg_9m",
                       "srad_seas_DJF","srad_seas_MAM","srad_seas_JJA","srad_seas_SON")
test_ext <- extent(c(clim_brick@extent@xmin+16000,
                     clim_brick@extent@xmin+17000,
                     clim_brick@extent@ymin+7000,
                     clim_brick@extent@ymin+11000))
test_env <- crop(clim_brick,test_ext,filename="~/GIS/microclim_preds_3m_test.grd",overwrite=TRUE)

##Loops through each model and generates seasonal raster predictions for disparity.
test_disp_stack <- stack()
for(i in 1:length(seas_disp_gam_mods)){
  print(paste("Now predicting raster for model",names(seas_disp_gam_mods)[i]))
  filename <- paste("~/GIS/test_",names(seas_disp_gam_mods)[i],".tif",sep="") 
  rast <- raster::predict(test_env,seas_disp_gam_mods[[i]],overwrite=TRUE,se.fit=TRUE,index=1:2,
                          filename=filename,const=data.frame(type="DataLogger"))
  names(rast) <- c(paste(names(seas_disp_gam_mods)[i],".est",sep=""),
                   paste(names(seas_disp_gam_mods)[i],".se",sep=""))
  test_disp_stack <- stack(test_disp_stack,rast)
}
plot(test_disp_stack)

test_sens_stack <- stack()
for(i in 1:length(seas_sens_gam_mods)){
  print(paste("Now predicting raster for model",names(seas_sens_gam_mods)[i]))
  filename <- paste("~/GIS/test_",names(seas_sens_gam_mods)[i],".tif",sep="") 
  rast <- raster::predict(test_env,seas_sens_gam_mods[[i]],overwrite=TRUE,se.fit=TRUE,index=1:2,
                          filename=filename,const=data.frame(type="DataLogger"))
  names(rast) <- c(paste(names(seas_sens_gam_mods)[i],".est",sep=""),
                   paste(names(seas_sens_gam_mods)[i],".se",sep=""))
  test_sens_stack <- stack(test_sens_stack,rast)
}
plot(test_sens_stack)

##Writes test stacks to disk.
test_disp <- brick(test_disp_stack,filename="~/GIS/test_disp.grd",overwrite=TRUE)
test_sens <- brick(test_sens_stack,filename="~/GIS/test_sens.grd",overwrite=TRUE)

##Loops through each model and generates seasonal raster predictions for disparity.
disp_stack <- stack()
for(i in 1:length(seas_disp_gam_mods)){
  print(paste("Now predicting raster for model",names(seas_disp_gam_mods)[i]))
  filename <- paste("/Volumes/ib_working/GIS/clim_",names(seas_disp_gam_mods)[i],".grd",sep="") 
  rast <- raster::predict(clim_brick,seas_disp_gam_mods[[i]],overwrite=TRUE,se.fit=TRUE,index=1:2,
                          filename=filename,datatype="FLT4S",const=data.frame(type="DataLogger"))
  names(rast) <- c(paste(names(seas_disp_gam_mods)[i],".est",sep=""),
                   paste(names(seas_disp_gam_mods)[i],".se",sep=""))
  disp_stack <- stack(disp_stack,rast)
}

sens_stack <- stack()
for(i in 1:length(seas_sens_gam_mods)){
  print(paste("Now predicting raster for model",names(seas_sens_gam_mods)[i]))
  filename <- paste("/Volumes/ib_working/GIS/clim_",names(seas_sens_gam_mods)[i],".grd",sep="") 
  rast <- raster::predict(clim_brick,seas_sens_gam_mods[[i]],overwrite=TRUE,se.fit=TRUE,index=1:2,
                          filename=filename,datatype="FLT4S",const=data.frame(type="DataLogger"))
  names(rast) <- c(paste(names(seas_sens_gam_mods)[i],".est",sep=""),
                   paste(names(seas_sens_gam_mods)[i],".se",sep=""))
  sens_stack <- stack(sens_stack,rast)
}

####Models of wet and dry season variables for paper####

##Scales UTM coordinates
meta_indices_tmax$utmx_scale <- as.numeric(scale(meta_indices_tmax$utmx))
meta_indices_tmax$utmy_scale <- as.numeric(scale(meta_indices_tmax$utmy))
meta_indices_tmin$utmx_scale <- as.numeric(scale(meta_indices_tmin$utmx))
meta_indices_tmin$utmy_scale <- as.numeric(scale(meta_indices_tmin$utmy))
meta_indices_tavg$utmx_scale <- as.numeric(scale(meta_indices_tavg$utmx))
meta_indices_tavg$utmy_scale <- as.numeric(scale(meta_indices_tavg$utmy))

##Fits the models for sensitivity
logit <- function(x) { exp(x) / (1-exp(x))}
tmax_gam <- gam(sens~ccov_81m*elev+I(srad_sum_9m/1e6)+utmx_scale+utmy_scale,data=meta_indices_tmax)
summary(tmax_gam)
#plot(tmax_gam)
tmin_gam <- gam(sens~s(elev,cold_ind_9m,k=5)+utmx_scale+utmy_scale,data=meta_indices_tmin)
summary(tmin_gam)
#plot(tmin_gam)
tavg_gam <- gam(sens~s(elev,cold_ind_9m,k=5)+I(sqrt(cvol_81m))+utmx_scale+utmy_scale,data=meta_indices_tavg)
summary(tavg_gam)
#plot(tavg_gam)

##Fits the models for disparity
disp_tmax_gam <- gam(disp~s(elev,cold_ind_9m,k=5)+s(cvol_81m,k=10)+I(srad_sum_9m/1e6)+utmx_scale+utmy_scale,data=meta_indices_tmax)
summary(disp_tmax_gam)
#plot(disp_tmax_gam)
disp_tmin_gam <- gam(disp~s(elev,cold_ind_9m,k=5)+utmx_scale+utmy_scale,data=meta_indices_tmin)
summary(disp_tmin_gam)
#plot(disp_tmin_gam)
disp_tavg_gam <- gam(disp~s(elev,cold_ind_9m,k=5)+ccov_81m+I(srad_sum_9m/1e6)+utmx_scale+utmy_scale,data=meta_indices_tavg)
summary(disp_tavg_gam)
#plot(disp_tavg_gam)

##Creates GAM predictions.
cold_ind_med <- median(meta_indices_tavg$cold_ind_9m)
can_vol_med <- median(meta_indices_tavg$cvol_81m)
srad_med <- median(meta_indices_tavg$srad_sum_9m)

pred_data_cold <- expand.grid(type="Cold-Air Index",
                              cold_ind_9m=seq(0.3,0.8,by=0.02),
                              elev=seq(500,2000,by=50),
                              cvol_81m=can_vol_med,
                              srad_sum_9m=srad_med,
                              utmx_scale=0,
                              utmy_scale=0)
pred_data_can <- expand.grid(type="Canopy Volume",
                             cold_ind_9m=cold_ind_med,
                             elev=seq(500,2000,by=50),
                             cvol_81m=seq(0,22000,by=1000),
                             srad_sum_9m=srad_med,
                             utmx_scale=0,
                             utmy_scale=0)
pred_data <- rbind(pred_data_cold,pred_data_can)

pred_data$Tmax_sen_gampred <- predict(tmax_gam,newdata=pred_data)
pred_data$Tmax_sen_gampred_se <- predict(tmax_gam,newdata=pred_data,se.fit=TRUE)$se.fit
pred_data$Tmin_sen_gampred <- predict(tmin_gam,newdata=pred_data)
pred_data$Tmin_sen_gampred_se <- predict(tmin_gam,newdata=pred_data,se.fit=TRUE)$se.fit
pred_data$Tavg_sen_gampred <- predict(tavg_gam,newdata=pred_data)
pred_data$Tavg_sen_gampred_se <- predict(tavg_gam,newdata=pred_data,se.fit=TRUE)$se.fit

##Estimates if sensitivities are significantly different from 1.
pred_data$Tmax_sen_gampred_low <- pred_data$Tmax_sen_gampred - pred_data$Tmax_sen_gampred_se * 2
pred_data$Tmax_sen_gampred_high <- pred_data$Tmax_sen_gampred + pred_data$Tmax_sen_gampred_se * 2
pred_data$Tmax_sen_gampred_sig <- as.numeric(!(pred_data$Tmax_sen_gampred_low < 1 & pred_data$Tmax_sen_gampred_high > 1))
pred_data$Tmin_sen_gampred_low <- pred_data$Tmin_sen_gampred - pred_data$Tmin_sen_gampred_se * 2
pred_data$Tmin_sen_gampred_high <- pred_data$Tmin_sen_gampred + pred_data$Tmin_sen_gampred_se * 2
pred_data$Tmin_sen_gampred_sig <- as.numeric(!(pred_data$Tmin_sen_gampred_low < 1 & pred_data$Tmin_sen_gampred_high > 1))
pred_data$Tavg_sen_gampred_low <- pred_data$Tavg_sen_gampred - pred_data$Tavg_sen_gampred_se * 2
pred_data$Tavg_sen_gampred_high <- pred_data$Tavg_sen_gampred + pred_data$Tavg_sen_gampred_se * 2
pred_data$Tavg_sen_gampred_sig <- as.numeric(!(pred_data$Tavg_sen_gampred_low < 1 & pred_data$Tavg_sen_gampred_high > 1))

##Plots predictions and standard errors.
p1 <- ggplot(subset(pred_data,type=="Cold-Air Index"))+
  geom_raster(aes(y=elev/1000,x=cold_ind_9m,fill=Tmax_sen_gampred))+
  geom_point(aes(y=elev/1000,
                 x=cold_ind_9m,
                 shape=factor(Tmax_sen_gampred_sig),
                 alpha=factor(Tmax_sen_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  scale_fill_gradient2(limits=c(0.8,1.2),midpoint=1,low="blue",high="red")+
  ylab("Elevation (km)")+
  xlab("Cold-Air Index (Unitless)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  #facet_wrap(facets=~type)+
  guides(fill = "none",shape="none",alpha="none")+
  theme_bw()
p2 <- ggplot(subset(pred_data,type=="Canopy Volume"))+
  geom_raster(aes(y=elev/1000,x=cvol_81m,fill=Tmax_sen_gampred))+
  scale_fill_gradient2("Climate \nExposure",limits=c(0.7,1.2),midpoint=1,low="blue",high="red")+
  geom_point(aes(y=elev/1000,
                 x=cvol_81m,
                 shape=factor(Tmax_sen_gampred_sig),
                 alpha=factor(Tmax_sen_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  ylab("")+
  xlab("Canopy Volume (m^3)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  #facet_wrap(facets=~type)+
  guides(shape="none",alpha="none")+
  theme_bw()+
  theme(axis.text.y = element_blank())
p3 <- ggplot(subset(pred_data,type=="Cold-Air Index"))+
  geom_raster(aes(y=elev/1000,x=cold_ind_9m,fill=Tmin_sen_gampred))+
  scale_fill_gradient2(limits=c(0.7,1.2),midpoint=1,low="blue",high="red")+
  geom_point(aes(y=elev/1000,
                 x=cold_ind_9m,
                 shape=factor(Tmin_sen_gampred_sig),
                 alpha=factor(Tmin_sen_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  ylab("Elevation (km)")+
  xlab("Cold-Air Index (Unitless)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  #facet_wrap(facets=~type)+
  guides(fill = "none",shape="none",alpha="none")+
  theme_bw()
p4 <- ggplot(subset(pred_data,type=="Canopy Volume"))+
  geom_raster(aes(y=elev/1000,x=cvol_81m,fill=Tmin_sen_gampred))+
  scale_fill_gradient2("Climate \nExposure",limits=c(0.7,1.2),midpoint=1,low="blue",high="red")+
  geom_point(aes(y=elev/1000,
                 x=cvol_81m,
                 shape=factor(Tmin_sen_gampred_sig),
                 alpha=factor(Tmin_sen_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  ylab("")+
  xlab("Canopy Volume (m^3)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  #facet_wrap(facets=~type)+
  guides(shape="none",alpha="none")+
  theme_bw()+
  theme(axis.text.y = element_blank())
p5 <- ggplot(subset(pred_data,type=="Cold-Air Index"))+
  geom_raster(aes(y=elev/1000,x=cold_ind_9m,fill=Tavg_sen_gampred))+
  scale_fill_gradient2(limits=c(0.7,1.2),midpoint=1,low="blue",high="red")+
  geom_point(aes(y=elev/1000,
                 x=cold_ind_9m,
                 shape=factor(Tavg_sen_gampred_sig),
                 alpha=factor(Tavg_sen_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  ylab("Elevation (km)")+
  xlab("Cold-Air Index (Unitless)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = "none",shape="none",alpha="none")+
  theme_bw()
p6 <- ggplot(subset(pred_data,type=="Canopy Volume"))+
  geom_raster(aes(y=elev/1000,x=cvol_81m,fill=Tavg_sen_gampred))+
  scale_fill_gradient2("Climate \nExposure",limits=c(0.7,1.2),midpoint=1,low="blue",high="red")+
  geom_point(aes(y=elev/1000,
                 x=cvol_81m,
                 shape=factor(Tavg_sen_gampred_sig),
                 alpha=factor(Tavg_sen_gampred_sig)),
             size=0.5)+
  scale_shape_manual(values=c(1,3))+
  scale_alpha_discrete(range=c(0,0.8))+
  ylab("")+
  xlab("Canopy Volume (m^3)")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  guides(shape = "none",alpha="none")+
  theme_bw()+
  theme(axis.text.y = element_blank())


pdf("../results/sum_tmax_sens_gam_2d.pdf",width=8,height=4)
grid.arrange(p1, p2, ncol=2, nrow=1,widths=c(3.5,4.5))
dev.off()

pdf("../results/win_tmin_sens_gam_2d.pdf",width=8,height=4)
grid.arrange(p3, p4, ncol=2, nrow=1,widths=c(3.5,4.5))
dev.off()

pdf("../results/tavg_sens_gam_2d.pdf",width=8,height=4)
grid.arrange(p5, p6, ncol=2, nrow=1,widths=c(3.5,4.5))
dev.off()

##Creates raster predictions.
library(raster)

##Prepares Rasters for prediction.
elev <- raster("/Volumes/ib_working/GIS/MORA_elev_3m.tif")
srad_sum_9m <- raster("/Volumes/ib_working/GIS/MORA_srad_yearsum_9m.tif")
utmx <- raster("/Volumes/ib_working/GIS/MORA_UTMX_3m.tif")
utmy <- raster("/Volumes/ib_working/GIS/MORA_UTMy_3m.tif")
cvol_81m <- raster("/Volumes/ib_working/GIS/MORA_can_vol_81m.tif")
ccov_81m <- raster("/Volumes/ib_working/GIS/MORA_can_pct_focal81m.tif")
cold_ind_9m <- raster("/Volumes/ib_working/GIS/MORA_coldair_3m.tif")
dry_ind_81m <- raster("/Volumes/ib_working/GIS/MORA_dry_index_81m_precip.tif")
clim_preds <- stack(elev,srad_sum_9m,utmx,utmy,cvol_81m,ccov_81m,cold_ind_9m,dry_ind_81m)

##Scaling for UTM coordinates.
utmx_center <- attributes(scale(meta_indices_tmax$utmx))$"scaled:center"
utmx_scale <- attributes(scale(meta_indices_tmax$utmx))$"scaled:scale"
utmy_center <- attributes(scale(meta_indices_tmax$utmy))$"scaled:center"
utmy_scale <- attributes(scale(meta_indices_tmax$utmy))$"scaled:scale"

clim_preds$utmx_scale <- (utmx - utmx_center) / utmx_scale
clim_preds$utmy_scale <- (utmy - utmy_center) / utmy_scale

clim_brick <- brick(clim_preds,format="raster",overwrite=TRUE,
                    filename="~/GIS/microclim_preds_3m.grd")

##Reads covariates.
#clim_brick <- brick("~/GIS/microclim_preds_3m.grd")
names(clim_brick) <- c("elev","srad_sum_9m","utmx","utmy","cvol_81m","ccov_81m",
                       "cold_ind_9m","dry_ind_81m","utmx_scale","utmy_scale")

##Subsets for testing.
##Tests the model on a small extent.
test_ext <- extent(c(clim_brick@extent@xmin+15000,
                     clim_brick@extent@xmin+19000,
                     clim_brick@extent@ymin+6000,
                     clim_brick@extent@ymin+10000))
test_env <- crop(clim_brick,test_ext,filename="~/GIS/microclim_preds_3m_test.grd")

##Predictions for test extent.
test_tmax_sens <- raster::predict(test_env,tmax_gam,overwrite=TRUE,
                                  filename="~/GIS/test_tmax_sens.tif")
test_tmax_disp <- raster::predict(test_env,disp_tmax_gam,overwrite=TRUE,
                                  filename="~/GIS/test_tmax_disp.tif")
test_tmin_sens <- raster::predict(test_env,tmin_gam,overwrite=TRUE,
                                  filename="~/GIS/test_tmin_sens.tif")
test_tmin_disp <- raster::predict(test_env,disp_tmin_gam,overwrite=TRUE,
                                  filename="~/GIS/test_tmin_disp.tif")
test_tavg_sens <- raster::predict(test_env,tavg_gam,overwrite=TRUE,
                                  filename="~/GIS/test_tavg_sens.tif")
test_tavg_disp <- raster::predict(test_env,disp_tavg_gam,overwrite=TRUE,
                                  filename="~/GIS/test_tavg_disp.tif")
par(mfrow=c(3,2))
plot(test_tmax_sens)
plot(test_tmax_disp)
plot(test_tmin_sens)
plot(test_tmin_disp)
plot(test_tavg_sens)
plot(test_tavg_disp)

##Predictions for full extent.
tmax_sens <- raster::predict(clim_brick,tmax_gam,overwrite=TRUE,progress="text",
                                  filename="~/GIS/tmax_sens_gam_pred.tif")
tmax_disp <- raster::predict(clim_brick,disp_tmax_gam,overwrite=TRUE,progress="text",
                                  filename="~/GIS/tmax_disp_gam_pred.tif")
tmin_sens <- raster::predict(clim_brick,tmin_gam,overwrite=TRUE,progress="text",
                                  filename="~/GIS/tmin_sens_gam_pred.tif")
tmin_disp <- raster::predict(clim_brick,disp_tmin_gam,overwrite=TRUE,progress="text",
                                  filename="~/GIS/tmin_disp_gam_pred.tif")
tavg_sens <- raster::predict(clim_brick,tavg_gam,overwrite=TRUE,progress="text",
                                  filename="~/GIS/tavg_sens_gam_pred.tif")
tavg_disp <- raster::predict(clim_brick,disp_tavg_gam,overwrite=TRUE,progress="text",
                                  filename="~/GIS/tavg_disp_gam_pred.tif")

##Standard errors for full extent.
tmax_sens_se <- raster::predict(clim_brick,tmax_gam,overwrite=TRUE,progress="text",
                             filename="~/GIS/tmax_sens_gam_se.tif",se.fit=TRUE,index=2)
tmax_disp_se <- raster::predict(clim_brick,disp_tmax_gam,overwrite=TRUE,progress="text",
                             filename="~/GIS/tmax_disp_gam_se.tif",se.fit=TRUE,index=2)
tmin_sens_se <- raster::predict(clim_brick,tmin_gam,overwrite=TRUE,progress="text",
                             filename="~/GIS/tmin_sens_gam_se.tif",se.fit=TRUE,index=2)
tmin_disp_se <- raster::predict(clim_brick,disp_tmin_gam,overwrite=TRUE,progress="text",
                             filename="~/GIS/tmin_disp_gam_se.tif",se.fit=TRUE,index=2)
tavg_sens_se <- raster::predict(clim_brick,tavg_gam,overwrite=TRUE,progress="text",
                             filename="~/GIS/tavg_sens_gam_se.tif",se.fit=TRUE,index=2)
tavg_disp_se <- raster::predict(clim_brick,disp_tavg_gam,overwrite=TRUE,progress="text",
                             filename="~/GIS/tavg_disp_gam_se.tif",se.fit=TRUE,index=2)
##PRISM regional means
prism_tavg <- 7.109838
prism_win_tmin <- -0.3824185
prism_sum_tmax <- 20.30168


