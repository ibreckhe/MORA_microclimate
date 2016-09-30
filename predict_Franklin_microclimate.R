##Script to generate predictions of microclimate variables at Franklin vegetation plot sites.

##Sets up workspace.
library(mgcv)
library(ggplot2)
library(dplyr)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")

##Reads in data.
tavg_meas <- read.csv("../results/buff_sens_decoup_tavg_covars.csv")
tmax_meas <- read.csv("../results/buff_sens_decoup_sum_tmax_covars.csv")
tmin_meas <- read.csv("../results/buff_sens_decoup_win_tmin_covars.csv")

fsites_covar <- read.csv("../results/franklin_ssmod_tavg_sens.csv")

prism <- read.csv("../raw/PRISM_daily_1_2004_9_2015.csv")
prism$DATE <- as.POSIXct(prism$DATE)
prism$TAVG <- (prism$TMAX + prism$TMIN) / 2

##Extracts PRISM regional estimates for the time period.
#Temporal bounds of the analysis
start_date <- as.POSIXct("2011-9-1")
end_date <- as.POSIXct("2015-9-1")

start_num_prism <- which(prism$DATE == start_date) + 1
end_num_prism <- which(prism$DATE == end_date) + 1

prism <- prism[start_num_prism:end_num_prism,]

##Separates winter and summer PRISM data.
prism$MONTH <- as.numeric(format(prism$DATE,format="%m"))
prism$SEASON <- NA
prism$SEASON[prism$MONTH %in% c(10,11,12,1,2,3,4,5)] <- "Winter"
prism$SEASON[prism$MONTH %in% c(6,7,8,9)] <- "Summer"

prism_sum <- prism
prism_sum[prism_sum$SEASON=="Winter",c(3:6)] <- NA

prism_win <- prism
prism_win[prism_win$SEASON=="Summer",c(3:6)] <- NA

##Computes yearly and seasonal averages.
prism_tavg <- mean(prism$TAVG)
prism_sum_tmax <- mean(prism_sum$TMAX,na.rm=TRUE)
prism_win_tmin <- mean(prism_win$TMIN,na.rm=TRUE)


##Scales UTM coordinates
tmax_meas$utmx_scale <- tmax_meas$utmx/1000
tmax_meas$utmy_scale <- tmax_meas$utmy/1000
tmin_meas$utmx_scale <- tmin_meas$utmx/1000
tmin_meas$utmy_scale <- tmin_meas$utmy/1000
tavg_meas$utmx_scale <- tavg_meas$utmx/1000
tavg_meas$utmy_scale <- tavg_meas$utmy/1000

##Filters out measurements with large bootstrap errors.
tmax_filt <- filter(tmax_meas,(sens_upr95 - sens_lwr95) < 0.1)
tmin_filt <- filter(tmin_meas,(sens_upr95 - sens_lwr95) < 0.1)
tavg_filt <- filter(tavg_meas,(sens_upr95 - sens_lwr95) < 0.05)

##Fits the models for sensitivity
tmax_gam <- gam(sens~s(elev,cold_ind_9m,k=10)+s(cvol_81m,k=5)+srad_sum_9m+utmx_scale+utmy_scale,data=tmax_filt)
summary(tmax_gam)
#plot(tmax_gam)
tmin_gam <- gam(sens~s(elev,cold_ind_9m,k=5)+utmx_scale+utmy_scale,data=tmin_filt)
summary(tmin_gam)
#plot(tmin_gam)
tavg_gam <- gam(sens~s(elev,cold_ind_9m,k=5)+cvol_81m+srad_sum_9m+utmx_scale+utmy_scale,data=tavg_filt)
summary(tavg_gam)
#plot(tavg_gam)

##Fits the models for disparity
disp_tmax_gam <- gam(disp~s(elev,cold_ind_9m,k=5)+s(cvol_81m,k=10)+srad_sum_9m+utmx_scale+utmy_scale,data=tmax_filt)
summary(disp_tmax_gam)
#plot(disp_tmax_gam)
disp_tmin_gam <- gam(disp~s(elev,cold_ind_9m,k=5)+utmx_scale+utmy_scale,data=tmin_filt)
summary(disp_tmin_gam)
#plot(disp_tmin_gam)
disp_tavg_gam <- gam(disp~s(elev,cold_ind_9m,k=5)+dry_ind_81m+utmx_scale+utmy_scale,data=tavg_filt)
summary(disp_tavg_gam)
#plot(disp_tavg_gam)

##Predicts all values at the Franklin locations.
fsites_pred <- data.frame(site=fsites_covar$SITE,
                          utmx_scale=fsites_covar$utmx/1000,
                          utmy_scale=fsites_covar$utmy/1000,
                          elev=fsites_covar$MORA_elev_3m,
                          cold_ind_9m=fsites_covar$MORA_coldair_index,
                          dry_ind_81m=fsites_covar$MORA_dry_index_81m_precip,
                          cvol_81m=fsites_covar$MORA_can_vol_81m,
                          srad_sum_9m=fsites_covar$MORA_srad_yearsum_9m)

fsites_pred$tmax_sens_est <- predict(tmax_gam,newdata=fsites_pred)
fsites_pred$tmax_sens_se <- predict(tmax_gam,newdata=fsites_pred,se.fit=TRUE)$se
fsites_pred$tmin_sens_est <- predict(tmin_gam,newdata=fsites_pred)
fsites_pred$tmin_sens_se <- predict(tmin_gam,newdata=fsites_pred,se.fit=TRUE)$se
fsites_pred$tavg_sens_est <- predict(tavg_gam,newdata=fsites_pred)
fsites_pred$tavg_sens_se <- predict(tavg_gam,newdata=fsites_pred,se.fit=TRUE)$se

fsites_pred$tmax_disp_est <- predict(disp_tmax_gam,newdata=fsites_pred)
fsites_pred$tmax_disp_se <- predict(disp_tmax_gam,newdata=fsites_pred,se.fit=TRUE)$se
fsites_pred$tmin_disp_est <- predict(disp_tmin_gam,newdata=fsites_pred)
fsites_pred$tmin_disp_se <- predict(disp_tmin_gam,newdata=fsites_pred,se.fit=TRUE)$se
fsites_pred$tavg_disp_est <- predict(disp_tavg_gam,newdata=fsites_pred)
fsites_pred$tavg_disp_se <- predict(disp_tavg_gam,newdata=fsites_pred,se.fit=TRUE)$se

##Joins with original covariates.
fsites_micro <- merge(fsites_covar,fsites_pred[,c(1,9:20)],by.x="SITE",by.y="site",keep.x=TRUE)

fsites_micro$tavg_est <- fsites_micro$tavg_disp_est + prism_tavg
fsites_micro$tmax_est <- fsites_micro$tmax_disp_est + prism_sum_tmax
fsites_micro$tmin_est <- fsites_micro$tmin_disp_est + prism_win_tmin

##Writes microclimate estimates to disk.
write.csv(fsites_micro,"../results/franklin_sites_microclimate.csv")
