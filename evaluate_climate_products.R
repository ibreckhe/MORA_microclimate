##Script to evaluate MORAclim products.

library(dismo)
library(raster)
library(dplyr)
source("~/code/MORA_microclimate/air_temp_functions.R")

####Loads daily temp data and site metadata.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
meta <- read.csv("/Users/ian/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/airtemp_sensor_metadata_2016.csv")
meta$alt_code <- gsub(pattern="-",replacement=".",x=meta$combined_1,fixed=TRUE)

prism <- read.csv("PRISM_daily_1_2004_9_2015.csv")
prism$DATE <- as.Date(prism$DATE)
days <- seq(as.Date("2009-09-01"),as.Date("2015-09-30"),by="day")
prism <- prism[prism$DATE %in% days,]
prism$DOY <- as.numeric(strftime(prism$DATE,format="%j"))
prism$MONTH <- as.numeric(strftime(prism$DATE,format="%m"))
prism$SEAS <- NA
prism$SEAS[prism$MONTH %in% c(12,1,2)] <- "DJF"
prism$SEAS[prism$MONTH %in% c(3,4,5)] <- "MAM"
prism$SEAS[prism$MONTH %in% c(6,7,8)] <- "JJA"
prism$SEAS[prism$MONTH %in% c(9,10,11)] <- "SON"

##Computes seasonal averages.
prism_gr <- group_by(prism,SEAS)
prism_seas <- summarise(prism_gr,avg_tmin=mean(TMIN),
                                 avg_tmax=mean(TMAX),
                                 avg_tavg=mean(TEMP),
                                 avg_prcp=mean(PRCP))


setwd("../cleaned/")
tavg_all <- read.csv("tavg_daily_all_stations_2009_2015.csv")
tavg_all$DATE <- as.Date(tavg_all$DATE)
tavg_all <- tavg_all[tavg_all$DATE %in% days,]
tmax_all <- read.csv("tmax_daily_all_stations_2009_2015.csv")
tmax_all$DATE <- as.Date(tmax_all$DATE)
tmax_all <- tmax_all[tmax_all$DATE %in% days,]
tmin_all <- read.csv("tmin_daily_all_stations_2009_2015.csv")
tmin_all$DATE <- as.Date(tmin_all$DATE)
tmin_all <- tmin_all[tmin_all$DATE %in% days,]

##Drops data from stations used for training.
stations_test <- c("AB08.STR.A1","AE10.A1","AM16.A1","AV14.A1","AX15.A2",
                   "HIGH.A1","LP19","OHAN.STR.A1","PA.1570m.A1","SPRY.STR.A2",
                   "SUNR.STR.A1","TB13.STR.A2","TO04.A1","TO04.A2",
                   "UpperPalisades","VATR.STR.A2","VATR.UPL.A2","WHRI.UPL.A2",
                   "WON8.A1","WON9.A1")

tavg_test <- tavg_all[,c("DATE",stations_test)]
tmin_test <- tmin_all[,c("DATE",stations_test)]
tmax_test <- tmax_all[,c("DATE",stations_test)]

##Subsets metadata at test locations
meta_test <- meta[meta$alt_code %in% stations_test,]

##Brings in prediction rasters.
disp_stack_tavg_90m <- brick("~/GIS/clim_seas_disp_tavg_gam_90m.grd")
disp_stack_tmin_90m <- brick("~/GIS/clim_seas_disp_tmin_gam_90m.grd")
disp_stack_tmax_90m <- brick("~/GIS/clim_seas_disp_tmax_gam_90m.grd")
sens_stack_tavg_90m <- brick("~/GIS/clim_seas_sens_tavg_gam_90m.grd")
sens_stack_tmin_90m <- brick("~/GIS/clim_seas_sens_tmin_gam_90m.grd")
sens_stack_tmax_90m <- brick("~/GIS/clim_seas_sens_tmax_gam_90m.grd")

disp_stack_90m <- stack(disp_stack_tavg_90m,disp_stack_tmax_90m,disp_stack_tmin_90m)
sens_stack_90m <- stack(sens_stack_tavg_90m,sens_stack_tmax_90m,sens_stack_tmin_90m)


disp_test_vals <- extract(disp_stack_90m,cbind(meta_test$X_UTM,meta_test$Y_UTM),
                          method='simple')
sens_test_vals <- extract(sens_stack_90m,cbind(meta_test$X_UTM,meta_test$Y_UTM),
                          method='bilinear')
test_vals <- cbind(disp_test_vals,sens_test_vals)
meta_test_vals <- data.frame(meta_test,test_vals)

##Computes predicted stats for each test station.
prism_seas_stats <- left_join(prism,prism_seas,by="SEAS")
prism_seas_stats$TMAX_cent <- prism_seas_stats$TMAX - prism_seas_stats$avg_tmax
prism_seas_stats$TMIN_cent <- prism_seas_stats$TMIN - prism_seas_stats$avg_tmin
prism_seas_stats$TAVG_cent <- prism_seas_stats$TEMP - prism_seas_stats$avg_tavg

##Loops through each test station and computes performance stats.
test_sum_all <- data.frame()
for (i in 1:nrow(meta_test_vals)){
  seas <- data.frame(SEAS=prism_seas_stats$SEAS)
  sens_DJF_tavg <- meta_test_vals[i,"sens_tavg_DJF"]
  sens_MAM_tavg <- meta_test_vals[i,"sens_tavg_MAM"]
  sens_JJA_tavg <- meta_test_vals[i,"sens_tavg_JJA"]
  sens_SON_tavg <- meta_test_vals[i,"sens_tavg_SON"]
  sens_seas_tavg <- c(sens_DJF_tavg,sens_MAM_tavg,sens_JJA_tavg,sens_SON_tavg)
  
  sens_DJF_tmin <- meta_test_vals[i,"sens_tmin_DJF"]
  sens_MAM_tmin <- meta_test_vals[i,"sens_tmin_MAM"]
  sens_JJA_tmin <- meta_test_vals[i,"sens_tmin_JJA"]
  sens_SON_tmin <- meta_test_vals[i,"sens_tmin_SON"]
  sens_seas_tmin <- c(sens_DJF_tmin,sens_MAM_tmin,sens_JJA_tmin,sens_SON_tmin)

  sens_DJF_tmax <- meta_test_vals[i,"sens_tmax_DJF"]
  sens_MAM_tmax <- meta_test_vals[i,"sens_tmax_MAM"]
  sens_JJA_tmax <- meta_test_vals[i,"sens_tmax_JJA"]
  sens_SON_tmax <- meta_test_vals[i,"sens_tmax_SON"]
  sens_seas_tmax <- c(sens_DJF_tmax,sens_MAM_tmax,sens_JJA_tmax,sens_SON_tmax)
  
  disp_DJF_tavg <- meta_test_vals[i,"disp_tavg_DJF"]
  disp_MAM_tavg <- meta_test_vals[i,"disp_tavg_MAM"]
  disp_JJA_tavg <- meta_test_vals[i,"disp_tavg_JJA"]
  disp_SON_tavg <- meta_test_vals[i,"disp_tavg_SON"]
  disp_seas_tavg <- c(disp_DJF_tavg,disp_MAM_tavg,disp_JJA_tavg,disp_SON_tavg)
  
  disp_DJF_tmin <- meta_test_vals[i,"disp_tmin_DJF"]
  disp_MAM_tmin <- meta_test_vals[i,"disp_tmin_MAM"]
  disp_JJA_tmin <- meta_test_vals[i,"disp_tmin_JJA"]
  disp_SON_tmin <- meta_test_vals[i,"disp_tmin_SON"]
  disp_seas_tmin <- c(disp_DJF_tmin,disp_MAM_tmin,disp_JJA_tmin,disp_SON_tmin)
  
  disp_DJF_tmax <- meta_test_vals[i,"disp_tmax_DJF"]
  disp_MAM_tmax <- meta_test_vals[i,"disp_tmax_MAM"]
  disp_JJA_tmax <- meta_test_vals[i,"disp_tmax_JJA"]
  disp_SON_tmax <- meta_test_vals[i,"disp_tmax_SON"]
  disp_seas_tmax <- c(disp_DJF_tmax,disp_MAM_tmax,disp_JJA_tmax,disp_SON_tmax)
  
  seas_df <- data.frame(SEAS=c("DJF","MAM","JJA","SON"),
                             sens_tavg=sens_seas_tavg,
                             sens_tmin=sens_seas_tmin,
                             sens_tmax=sens_seas_tmax,
                             disp_tavg=disp_seas_tavg,
                             disp_tmin=disp_seas_tmin,
                             disp_tmax=disp_seas_tmax)
  prism_site_seas <- left_join(prism_seas_stats,seas_df,by="SEAS")
  prism_site_tavg <- (prism_site_seas$TAVG_cent * prism_site_seas$sens_tavg) + prism_site_seas$disp_tavg + prism_site_seas$avg_tavg
  prism_site_tmax <- (prism_site_seas$TMAX_cent * prism_site_seas$sens_tmax) + prism_site_seas$disp_tmax + prism_site_seas$avg_tmax
  prism_site_tmin <- (prism_site_seas$TMIN_cent * prism_site_seas$sens_tmin) + prism_site_seas$disp_tmin + prism_site_seas$avg_tmin
  obs_site_tavg <- tavg_test[,meta_test_vals$alt_code[i]]
  obs_site_tmin <- tmin_test[,meta_test_vals$alt_code[i]]
  obs_site_tmax <- tmax_test[,meta_test_vals$alt_code[i]]
  test_df <- data.frame(days,seas,prism_site_tavg,prism_site_tmax,prism_site_tmin,
                        obs_site_tavg,obs_site_tmin,obs_site_tmax)
  test_df <- test_df[complete.cases(test_df),]
  test_df_seas <- group_by(test_df,SEAS)
  test_sum_seas <- summarise(test_df_seas,
                    tmax_daily_mae=mean(abs(prism_site_tmax-obs_site_tmax)),
                    tmin_daily_mae=mean(abs(prism_site_tmin-obs_site_tmin)),
                    tavg_daily_mae=mean(abs(prism_site_tavg-obs_site_tavg)),
                    tmax_avg_err=mean(prism_site_tmax) - mean(obs_site_tmax),
                    tmin_avg_err=mean(prism_site_tmin) - mean(obs_site_tmin),
                    tavg_avg_err=mean(prism_site_tavg) - mean(obs_site_tavg),
                    tmax_obs_avg=mean(obs_site_tmax),
                    tmin_obs_avg=mean(obs_site_tmin),
                    tavg_obs_avg=mean(obs_site_tavg),
                    tmax_pred_avg=mean(prism_site_tmax),
                    tmin_pred_avg=mean(prism_site_tmin),
                    tavg_pred_avg=mean(prism_site_tavg))
  test_df_ann <- group_by(test_df)
  test_sum_ann <- summarise(test_df_ann,
                             tmax_daily_mae=mean(abs(prism_site_tmax-obs_site_tmax)),
                             tmin_daily_mae=mean(abs(prism_site_tmin-obs_site_tmin)),
                             tavg_daily_mae=mean(abs(prism_site_tavg-obs_site_tavg)),
                             tmax_avg_err=mean(prism_site_tmax) - mean(obs_site_tmax),
                             tmin_avg_err=mean(prism_site_tmin) - mean(obs_site_tmin),
                             tavg_avg_err=mean(prism_site_tavg) - mean(obs_site_tavg),
                             tmax_obs_avg=mean(obs_site_tmax),
                             tmin_obs_avg=mean(obs_site_tmin),
                             tavg_obs_avg=mean(obs_site_tavg),
                             tmax_pred_avg=mean(prism_site_tmax),
                             tmin_pred_avg=mean(prism_site_tmin),
                             tavg_pred_avg=mean(prism_site_tavg))
                             
  test_sum_ann$SEAS <- "Annual"
  test_sum <- rbind(test_sum_seas,test_sum_ann)
  test_sum$site <- meta_test_vals$alt_code[i]
  test_sum_all <- rbind(test_sum_all,test_sum)
}

###Plots performance stats.
library(ggplot2)
library(reshape2)
library(gridExtra)

plot(density(test_sum_all$tmin_avg_err))
plot(density(test_sum_all$tmax_avg_err))
plot(density(test_sum_all$tavg_avg_err))
test_sum_melt <- melt(test_sum_all)
test_sum_melt$variable <- factor(test_sum_melt$variable,
                                labels=c("Tmax MAE","Tmin MAE","Tavg MAE",
                                         "Tmax","Tmin","Tavg",
                                         "Obs. Tmax","Obs. Tmin","Obs. Tavg",
                                         "Pred. Tmax","Pred. Tmin","Pred. Tavg"))
test_sum_melt$SEAS <- factor(test_sum_melt$SEAS,levels=c("DJF","MAM","JJA","SON","Annual"))

##Error histograms.
pdf("../results/MORAclim_test_averages_seasonal.pdf",width=4,height=5)
ggplot(subset(test_sum_melt,variable %in% c("Tmax","Tmin","Tavg")))+
  geom_histogram(aes(x=value),fill="grey60",color="white",
                 breaks=seq(-2.25,2.25,by=0.5))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  scale_y_continuous(breaks=seq(0,12,by=4))+
  facet_grid(facets=SEAS~variable)+
  theme_bw()+
  theme(panel.grid=element_blank())
dev.off()

##Observed vs Predicted Plots.
p1 <- ggplot(test_sum_all)+
        geom_point(aes(x=tmax_pred_avg,y=tmax_obs_avg,color=SEAS,shape=SEAS))+
        geom_abline(aes(slope=1,intercept=0),linetype="dotted")+
        xlab("")+
        scale_color_grey("Season",guide=FALSE)+
        scale_shape_discrete("Season",guide=FALSE)+
        ylab("Predicted")+
        ggtitle("Tmax")+
        theme_bw()+
        theme(panel.grid=element_blank())
p2 <- ggplot(test_sum_all)+
  geom_point(aes(x=tmin_pred_avg,y=tmin_obs_avg,color=SEAS,shape=SEAS))+
  geom_abline(aes(slope=1,intercept=0),linetype="dotted")+
  xlab("Observed")+
  scale_color_grey("Season",guide=FALSE)+
  scale_shape_discrete("Season",guide=FALSE)+
  ylab("")+
  ggtitle("Tmin")+
  theme_bw()+
  theme(panel.grid=element_blank())

p3 <- ggplot(test_sum_all)+
  geom_point(aes(x=tavg_pred_avg,y=tavg_obs_avg,color=SEAS,shape=SEAS))+
  geom_abline(aes(slope=1,intercept=0),linetype="dotted")+
  xlab("")+
  scale_color_grey("Season")+
  scale_shape_discrete("Season")+
  ylab("")+
  ggtitle("Tavg")+
  theme_bw()+
  theme(panel.grid=element_blank())

pdf("../results/MORAclim_test_obs_pred.pdf",width=10,height=3.5)
grid.arrange(p1,p2,p3,ncol=3,widths=c(1,1,1.3))
dev.off()

####Tests Kriged estimates of daily data.
setwd("/Volumes/ib_working/MORAclim_krige/")
est <- list.files(".",pattern="est.tif$")
est_tmin <- grep(est,pattern="tmin",fixed=TRUE,value=TRUE)
est_tmax <- grep(est,pattern="tmax",fixed=TRUE,value=TRUE)
est_tavg <- grep(est,pattern="tavg",fixed=TRUE,value=TRUE)

days <- seq(as.Date("2009-01-01"),as.Date("2015-09-30"),by="day")

##Creates an empty array to hold output.
test <- array(NA,dim=c(length(meta_test$alt_code),11,length(which(days==as.Date("2009-09-01")):length(days))))

daystart <- which(days==as.Date("2009-09-01"))
dayend <- which(days==as.Date("2015-09-30"))
dayseq <- daystart:dayend
daylab <- as.character(seq(as.Date("2009-09-01"),as.Date("2015-09-30"),by="day"))

for(i in 1:length(daystart:dayend)){
  print(paste("Processing ",days[dayseq[i]]))
  est_tmin_day <- raster(est_tmin[dayseq[i]])
  est_tmax_day <- raster(est_tmax[dayseq[i]])
  est_tavg_day <- raster(est_tavg[dayseq[i]])
  
  obs_tmin_day <- t(tmin_test[tmin_test$DATE == days[dayseq[i]],][,-1])
  obs_tmax_day <- t(tmax_test[tmax_test$DATE == days[dayseq[i]],][,-1])
  obs_tavg_day <- t(tavg_test[tavg_test$DATE == days[dayseq[i]],][,-1])
  
  obs_df <- data.frame(site=rownames(obs_tmin_day),
                       obs_tmin=obs_tmin_day[,1],
                       obs_tmax=obs_tmax_day[,1],
                       obs_tavg=obs_tavg_day[,1])
  meta_test_geo <- data.frame(site=meta_test$alt_code,
                              utmx=meta_test$X_UTM,
                              utmy=meta_test$Y_UTM)
  obs_df_geo <- left_join(obs_df,meta_test_geo,by="site")
  obs_df_geo$est_tmin <- extract(est_tmin_day,cbind(obs_df_geo$utmx,obs_df_geo$utmy),
                                                    method="bilinear")/10000
  obs_df_geo$est_tmax <- extract(est_tmax_day,cbind(obs_df_geo$utmx,obs_df_geo$utmy),
                                 method="bilinear")/10000
  obs_df_geo$est_tavg <- extract(est_tavg_day,cbind(obs_df_geo$utmx,obs_df_geo$utmy),
                                 method="bilinear")/10000
  obs_df_geo$tmin_err <- obs_df_geo$obs_tmin - obs_df_geo$est_tmin
  obs_df_geo$tmax_err <- obs_df_geo$obs_tmax - obs_df_geo$est_tmax
  obs_df_geo$tavg_err <- obs_df_geo$obs_tavg - obs_df_geo$est_tavg
  test[,,i] <- as.matrix(obs_df_geo[,-1])
}
dimnames(test) <- list(obs_df_geo$site,
                       colnames(obs_df_geo)[-1],
                       daylab)

##Summarises daily errors by season.
site_sums <- data.frame()
for(i in 1:length(obs_df_geo$site)){
  site_data <- data.frame(t(test[i,,]))
  site_data$date <- as.Date(rownames(site_data))
  site_data$month <- as.numeric(strftime(site_data$date,format="%m"))
  site_data$seas <- NA
  site_data$seas[site_data$month %in% c(12,1,2)] <- "DJF"
  site_data$seas[site_data$month %in% c(3,4,5)] <- "MAM"
  site_data$seas[site_data$month %in% c(6,7,8)] <- "JJA"
  site_data$seas[site_data$month %in% c(9,10,11)] <- "SON"
  site_data <- site_data[complete.cases(site_data),]
  site_data <- group_by(site_data,seas)
  site_sum <- summarise(site_data,
                        obs_tmin_mean=mean(obs_tmin),
                        obs_tmax_mean=mean(obs_tmax),
                        obs_tavg_mean=mean(obs_tavg),
                        est_tmin_mean=mean(est_tmin),
                        est_tmax_mean=mean(est_tmax),
                        est_tavg_mean=mean(est_tavg),
                        abs_err_tmin=mean(abs(tmin_err)),
                        abs_err_tmax=mean(abs(tmax_err)),
                        abs_err_tavg=mean(abs(tavg_err)),
                        min_err_tmin=min(tmin_err),
                        min_err_tmax=min(tmax_err),
                        min_err_tavg=min(tavg_err),
                        max_err_tmin=max(tmin_err),
                        max_err_tmax=max(tmax_err),
                        max_err_tavg=max(tavg_err))
  site_sum$site <- obs_df_geo$site[i]
  site_sums <- rbind(site_sums,site_sum)
}

##Density plots of daily errors.
pdf("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/MORAclim_daily_krige_err.pdf",width=9,height=4)
par(mfrow=c(1,3),oma=c(2,2,0,0),mar=c(2,2,1,1))
plot(density(test[1,9,],na.rm=TRUE),xlim=c(-5,5),ylim=c(0,2),main="Tmin",
     xlab="")
for(i in 2:20){
  points(density(test[i,9,],na.rm=TRUE),type="l",col=i)
}
plot(density(test[1,10,],na.rm=TRUE),xlim=c(-5,5),ylim=c(0,2),main="Tmax",
     xlab="Error (C)")
for(i in 2:20){
  points(density(test[i,10,],na.rm=TRUE),type="l",col=i)
}
plot(density(test[1,11,],na.rm=TRUE),xlim=c(-5,5),ylim=c(0,2),main="Tavg",
     xlab="")
for(i in 2:20){
  points(density(test[i,11,],na.rm=TRUE),type="l",col=i)
}
dev.off()

##Identifies days with large errors across sites.
abs_err_fun <- function(x){
  mean(abs(x),na.rm=TRUE)
}
avg_err_tmin_days <- apply(test[,9,],MARGIN=2,FUN=abs_err_fun)
avg_err_tmax_days <- apply(test[,10,],MARGIN=2,FUN=abs_err_fun)
avg_err_tavg_days <- apply(test[,11,],MARGIN=2,FUN=abs_err_fun)

##Plots time-series of average errors.
par(mfrow=c(7,1),oma=c(1,1,0,0))
for(i in 2009:2015){
  plot(as.Date(daylab),avg_err_tmin_days,type="l",ylim=c(0,3),col=4,
       xlim=as.Date(paste(i,c("-01-01","-12-31"),sep="")),main=as.character(i))
  points(as.Date(daylab),avg_err_tmax_days,type="l",col=2)
  points(as.Date(daylab),avg_err_tavg_days,type="l",col=1)
}

##Finds days with large errors.
prob_days <- daylab[which(avg_err_tmin_days > 1.5 | avg_err_tmax_days > 1.5 | avg_err_tavg_days > 1.5)]
prob_df <- data.frame(date=daylab,
                      avg_err_tmin=avg_err_tmin_days,
                      avg_err_tmax=avg_err_tmax_days,
                      avg_err_tavg=avg_err_tavg_days)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/")
write.csv(prob_df,"./results/Krige_test_err_days.csv",row.names=FALSE)

##Evaluates Bayesian temperature predictions.

setwd("/Volumes/ib_working/MORAclim_bayes/")
est <- list.files(".",pattern="est.tif$")
est_tmin <- grep(est,pattern="tmin",fixed=TRUE,value=TRUE)
est_tmax <- grep(est,pattern="tmax",fixed=TRUE,value=TRUE)
est_tavg <- grep(est,pattern="tavg",fixed=TRUE,value=TRUE)

##Creates an empty array to hold output.
test_bayes <- array(NA,dim=c(length(meta_test$alt_code),11,length(which(days==as.Date("2009-09-01")):length(days))))

daystart <- which(days==as.Date("2009-09-01"))
dayend <- which(days==as.Date("2015-09-30"))
dayseq <- daystart:dayend
daylab <- as.character(seq(as.Date("2009-09-01"),as.Date("2015-09-30"),by="day"))

for(i in 1:length(daystart:dayend)){
  print(paste("Processing ",days[dayseq[i]]))
  est_tmin_day <- raster(est_tmin[dayseq[i]])
  est_tmax_day <- raster(est_tmax[dayseq[i]])
  est_tavg_day <- raster(est_tavg[dayseq[i]])
  
  obs_tmin_day <- t(tmin_test[tmin_test$DATE == days[dayseq[i]],][,-1])
  obs_tmax_day <- t(tmax_test[tmax_test$DATE == days[dayseq[i]],][,-1])
  obs_tavg_day <- t(tavg_test[tavg_test$DATE == days[dayseq[i]],][,-1])
  
  obs_df <- data.frame(site=rownames(obs_tmin_day),
                       obs_tmin=obs_tmin_day[,1],
                       obs_tmax=obs_tmax_day[,1],
                       obs_tavg=obs_tavg_day[,1])
  meta_test_geo <- data.frame(site=meta_test$alt_code,
                              utmx=meta_test$X_UTM,
                              utmy=meta_test$Y_UTM)
  obs_df_geo <- left_join(obs_df,meta_test_geo,by="site")
  obs_df_geo$est_tmin <- extract(est_tmin_day,cbind(obs_df_geo$utmx,obs_df_geo$utmy),
                                 method="bilinear")/10000
  obs_df_geo$est_tmax <- extract(est_tmax_day,cbind(obs_df_geo$utmx,obs_df_geo$utmy),
                                 method="bilinear")/10000
  obs_df_geo$est_tavg <- extract(est_tavg_day,cbind(obs_df_geo$utmx,obs_df_geo$utmy),
                                 method="bilinear")/10000
  obs_df_geo$tmin_err <- obs_df_geo$obs_tmin - obs_df_geo$est_tmin
  obs_df_geo$tmax_err <- obs_df_geo$obs_tmax - obs_df_geo$est_tmax
  obs_df_geo$tavg_err <- obs_df_geo$obs_tavg - obs_df_geo$est_tavg
  test_bayes[,,i] <- as.matrix(obs_df_geo[,-1])
}
dimnames(test_bayes) <- list(obs_df_geo$site,
                       colnames(obs_df_geo)[-1],
                       daylab)

##Summarises daily errors by season.
site_sums_bayes <- data.frame()
for(i in 1:length(obs_df_geo$site)){
  site_data <- data.frame(t(test_bayes[i,,]))
  site_data$date <- as.Date(rownames(site_data))
  site_data$month <- as.numeric(strftime(site_data$date,format="%m"))
  site_data$seas <- NA
  site_data$seas[site_data$month %in% c(12,1,2)] <- "DJF"
  site_data$seas[site_data$month %in% c(3,4,5)] <- "MAM"
  site_data$seas[site_data$month %in% c(6,7,8)] <- "JJA"
  site_data$seas[site_data$month %in% c(9,10,11)] <- "SON"
  site_data <- site_data[complete.cases(site_data),]
  site_data <- group_by(site_data,seas)
  site_sum <- summarise(site_data,
                        obs_tmin_mean=mean(obs_tmin),
                        obs_tmax_mean=mean(obs_tmax),
                        obs_tavg_mean=mean(obs_tavg),
                        est_tmin_mean=mean(est_tmin),
                        est_tmax_mean=mean(est_tmax),
                        est_tavg_mean=mean(est_tavg),
                        abs_err_tmin=mean(abs(tmin_err)),
                        abs_err_tmax=mean(abs(tmax_err)),
                        abs_err_tavg=mean(abs(tavg_err)),
                        min_err_tmin=min(tmin_err),
                        min_err_tmax=min(tmax_err),
                        min_err_tavg=min(tavg_err),
                        max_err_tmin=max(tmin_err),
                        max_err_tmax=max(tmax_err),
                        max_err_tavg=max(tavg_err))
  site_sum$site <- obs_df_geo$site[i]
  site_sums_bayes <- rbind(site_sums_bayes,site_sum)
}

##Density plots of daily errors.
pdf("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/MORAclim_daily_bayes_err.pdf",width=9,height=4)
par(mfrow=c(1,3),oma=c(2,2,0,0),mar=c(2,2,1,1))
plot(density(test_bayes[1,9,],na.rm=TRUE),xlim=c(-10,10),ylim=c(0,2),main="Tmin",
     xlab="")
for(i in 2:20){
  points(density(test_bayes[i,9,],na.rm=TRUE),type="l",col=i)
}
plot(density(test_bayes[1,10,],na.rm=TRUE),xlim=c(-10,10),ylim=c(0,2),main="Tmax",
     xlab="Error (C)")
for(i in 2:20){
  points(density(test_bayes[i,10,],na.rm=TRUE),type="l",col=i)
}
plot(density(test_bayes[1,11,],na.rm=TRUE),xlim=c(-10,10),ylim=c(0,2),main="Tavg",
     xlab="")
for(i in 2:20){
  points(density(test_bayes[i,11,],na.rm=TRUE),type="l",col=i)
}
dev.off()

##Plots mean absolute error across seasons.
sites_bayes_melt <- melt(site_sums_bayes)
sites_bayes_err <- subset(sites_bayes_melt,variable %in% c("abs_err_tmax","abs_err_tmin","abs_err_tavg"))
sites_bayes_err$variable <- factor(sites_bayes_err$variable,
                                   levels=c("abs_err_tmax","abs_err_tavg","abs_err_tmin"),
                                   labels=c("Tmax","Tavg","Tmin"))
sites_bayes_err$seas <- factor(sites_bayes_err$seas,
                                   levels=c("DJF","MAM","JJA","SON"))
pdf("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/MORAclim_bayes_abs_err_seas.pdf",
    width=4,height=5)
ggplot(sites_bayes_err)+
  geom_boxplot(aes(x=seas,y=value,fill=variable),
               outlier.size=0,width=0.5,position=position_dodge(0.5))+
  scale_fill_grey("",start=0.4,end=0.9)+
  ylim(c(0,2))+
  ylab("Mean Absolute Error (C)")+
  xlab("Season")+
  theme_bw()+
  theme(panel.grid=element_blank())
dev.off()


