##Script to test temporal consistency of macroclimate-microclimate relationships.

source("~/code/MORA_microclimate/air_temp_functions.R")
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/")
library(lme4)
library(ggplot2)
library(dplyr)

##Loads site metadata.
meta <- read.csv("./raw/airtemp_sensor_metadata_2016.csv")
meta$site <- gsub("-",".",meta$combined_1,fixed=TRUE)
meta$site <- gsub("1791.STR","X1791.STR",meta$site,fixed=TRUE)
meta$site <- gsub("OhanapecoshNPS","Ohanapecosh",meta$site,fixed=TRUE)

##Loads tavg_data
tavg <- read.csv("./cleaned/yearly_seas_indices_tavg_2009_2015.csv")
tavg_meta <- left_join(tavg,meta,by="site")
tavg_meta$disp_range <- tavg_meta$disp_upr95 - tavg_meta$disp_lwr95
tavg_meta$sens_range <- tavg_meta$sens_upr95 - tavg_meta$sens_lwr95
tavg_meta$srad_seas <- NA
tavg_meta$srad_seas[tavg_meta$seas=="DJF"] <- tavg_meta$srad_seas_DJF[tavg_meta$seas=="DJF"]
tavg_meta$srad_seas[tavg_meta$seas=="MAM"] <- tavg_meta$srad_seas_MAM[tavg_meta$seas=="MAM"]
tavg_meta$srad_seas[tavg_meta$seas=="JJA"] <- tavg_meta$srad_seas_JJA[tavg_meta$seas=="JJA"]
tavg_meta$srad_seas[tavg_meta$seas=="SON"] <- tavg_meta$srad_seas_SON[tavg_meta$seas=="SON"]

##Loads tmin_data
tmin <- read.csv("./cleaned/yearly_seas_indices_tmin_2009_2015.csv")
tmin_meta <- left_join(tmin,meta,by="site")
tmin_meta$disp_range <- tmin_meta$disp_upr95 - tmin_meta$disp_lwr95
tmin_meta$sens_range <- tmin_meta$sens_upr95 - tmin_meta$sens_lwr95
tmin_meta$srad_seas <- NA
tmin_meta$srad_seas[tmin_meta$seas=="DJF"] <- tmin_meta$srad_seas_DJF[tmin_meta$seas=="DJF"]
tmin_meta$srad_seas[tmin_meta$seas=="MAM"] <- tmin_meta$srad_seas_MAM[tmin_meta$seas=="MAM"]
tmin_meta$srad_seas[tmin_meta$seas=="JJA"] <- tmin_meta$srad_seas_JJA[tmin_meta$seas=="JJA"]
tmin_meta$srad_seas[tmin_meta$seas=="SON"] <- tmin_meta$srad_seas_SON[tmin_meta$seas=="SON"]

##Loads tmax data
tmax <- read.csv("./cleaned/yearly_seas_indices_tmax_2009_2015.csv")
tmax_meta <- left_join(tmax,meta,by="site")
tmax_meta$disp_range <- tmax_meta$disp_upr95 - tmax_meta$disp_lwr95
tmax_meta$sens_range <- tmax_meta$sens_upr95 - tmax_meta$sens_lwr95
tmax_meta$srad_seas <- NA
tmax_meta$srad_seas[tmax_meta$seas=="DJF"] <- tmax_meta$srad_seas_DJF[tmax_meta$seas=="DJF"]
tmax_meta$srad_seas[tmax_meta$seas=="MAM"] <- tmax_meta$srad_seas_MAM[tmax_meta$seas=="MAM"]
tmax_meta$srad_seas[tmax_meta$seas=="JJA"] <- tmax_meta$srad_seas_JJA[tmax_meta$seas=="JJA"]
tmax_meta$srad_seas[tmax_meta$seas=="SON"] <- tmax_meta$srad_seas_SON[tmax_meta$seas=="SON"]


## Filters data to remove suspect measurements.
disp_thresh <- 1.3
disp_min <- -8
disp_max <- 8
sens_thresh <- 0.3
sens_min <- 0.65
sens_max <- 1.45

tavg_meta_sub <- filter(tavg_meta,disp_range < disp_thresh & 
                          sens_range < sens_thresh &
                          sens > sens_min &
                          sens < sens_max &
                          disp > disp_min &
                          disp < disp_max)
tmin_meta_sub <- filter(tmin_meta,disp_range < disp_thresh & 
                          sens_range < sens_thresh &
                          sens > sens_min &
                          sens < sens_max &
                          disp > disp_min &
                          disp < disp_max)
tmax_meta_sub <- filter(tmax_meta,disp_range < disp_thresh & 
                          sens_range < sens_thresh &
                          sens > sens_min &
                          sens < sens_max &
                          disp > disp_min &
                          disp < disp_max)
##Exploratory Plots

##Tavg
pdf("./results/tavg_disparity_5yr.pdf",width=7,height=5)
ggplot(tavg_meta_sub)+
  geom_point(aes(x=elev,y=disp,col=seas),size=0.5)+
  xlab("Elevation (m)")+
  scale_color_discrete("Season")+
  ylab("Local Disparity (C)")+
  geom_smooth(aes(x=elev,y=disp,col=seas),method="gam")+
  facet_wrap(facets=~year,nrow=2)+
  theme_bw()
dev.off()
ggplot(tavg_meta_sub)+
  geom_point(aes(x=elev,y=sens,col=seas))+
  geom_smooth(aes(x=elev,y=sens,col=seas),method="gam")+
  facet_grid(facets=.~year)+
  theme_bw()

#Tmin
ggplot(tmin_meta_sub)+
  geom_point(aes(x=elev,y=disp,col=seas))+
  geom_smooth(aes(x=elev,y=disp,col=seas),method="gam")+
  facet_grid(facets=.~year)+
  theme_bw()
ggplot(tmin_meta_sub)+
  geom_point(aes(x=elev,y=sens,col=seas))+
  geom_smooth(aes(x=elev,y=sens,col=seas),method="gam")+
  facet_grid(facets=.~year)+
  theme_bw()

#Tmax
ggplot(tmax_meta_sub)+
  geom_point(aes(x=elev,y=disp,col=seas))+
  geom_smooth(aes(x=elev,y=disp,col=seas),method="gam")+
  facet_grid(facets=.~year)+
  theme_bw()
ggplot(tmax_meta_sub)+
  geom_point(aes(x=elev,y=sens,col=seas))+
  geom_smooth(aes(x=elev,y=sens,col=seas),method="gam")+
  facet_grid(facets=.~year)+
  theme_bw()


##Initial model of tavg
tavg_disp_lm1 <- lm(disp~seas+elev+cair_reg_9m+X_UTM+Y_UTM+I(reg_cancov_243m/100)+
                      I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                      elev:seas+type,
                    data=tavg_meta_sub)
summary(tavg_disp_lm1)

tavg_sens_lm1 <- lm(sens~seas+elev+cair_reg_9m+X_UTM+Y_UTM+I(reg_cancov_243m/100)+cair_reg_9m:seas+elev:seas+type,
                    data=tavg_meta_sub)
summary(tavg_sens_lm1)

tmin_disp_lm1 <- lm(disp~seas+I(elev/1000)+cair_reg_9m+I(logit(reg_cancov_81m/100+0.05))+I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type,
                    data=tmin_meta_sub)
summary(tmin_disp_lm1)

tmin_sens_lm1 <- lm(sens~seas+I(elev/1000)+cair_reg_9m+I(logit(reg_cancov_81m/100+0.05))+I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type,
                    data=tmin_meta_sub)
summary(tmin_sens_lm1)

tmax_disp_lm1 <- lm(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                      I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type,
                    data=tmax_meta_sub)
summary(tmax_disp_lm1)

tmax_sens_lm1 <- lm(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                      I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type,
                    data=tmax_meta_sub)
summary(tmax_sens_lm1)

##Initial mixed-effects models without year effects.
tavg_disp_lme1 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+(1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme1)

tavg_sens_lme1 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+(1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme1)

tmin_disp_lme1 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+I(reg_cancov_243m/100):elev+I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+(1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme1)

tmin_sens_lme1 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+(1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme1)

tmax_disp_lme1 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+(1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme1)

tmax_sens_lme1 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+(1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme1)

##Tests the effect of adding a year term
tavg_disp_lme2 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme2)
AIC(tavg_disp_lme1,tavg_disp_lme2)

tavg_sens_lme2 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme2)
AIC(tavg_sens_lme1,tavg_sens_lme2)

tmin_disp_lme2 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme2)
AIC(tmin_disp_lme1,tmin_disp_lme2)

tmin_sens_lme2 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme2)
AIC(tmin_sens_lme1,tmin_sens_lme2)

tmax_disp_lme2 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme2)
AIC(tmax_disp_lme1,tmax_disp_lme2)

tmax_sens_lme2 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+
                       (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme2)
AIC(tmax_sens_lme1,tmax_sens_lme2)

##Adds a year:season interaction.
tavg_disp_lme3 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme3)
AIC(tavg_disp_lme1,tavg_disp_lme2,tavg_disp_lme3)

tavg_sens_lme3 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme3)
AIC(tavg_sens_lme1,tavg_sens_lme2,tavg_sens_lme3)

tmin_disp_lme3 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme3)
AIC(tmin_disp_lme1,tmin_disp_lme2,tmin_disp_lme3)

tmin_sens_lme3 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)
                       +I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme3)
AIC(tmin_sens_lme1,tmin_sens_lme2,tmin_sens_lme3)

tmax_disp_lme3 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme3)
AIC(tmax_disp_lme1,tmax_disp_lme2,tmax_disp_lme3)

tmax_sens_lme3 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme3)
AIC(tmax_sens_lme1,tmax_sens_lme2,tmax_sens_lme3)

##Adds a year:elev interaction.
tavg_disp_lme4 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+year:I(elev/1000)+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme4)
AIC(tavg_disp_lme1,tavg_disp_lme2,tavg_disp_lme3,tavg_disp_lme4)

tavg_sens_lme4 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+year:I(elev/1000)+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme4)
AIC(tavg_sens_lme1,tavg_sens_lme2,tavg_sens_lme3,tavg_sens_lme4)

tmin_disp_lme4 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+year:I(elev/1000)+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme4)
AIC(tmin_disp_lme1,tmin_disp_lme2,tmin_disp_lme3,tmin_disp_lme4)

tmin_sens_lme4 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)
                       +I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+year:I(elev/1000)+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme4)
AIC(tmin_sens_lme1,tmin_sens_lme2,tmin_sens_lme3,tmin_sens_lme4)

tmax_disp_lme4 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+year:I(elev/1000)+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme4)
AIC(tmax_disp_lme1,tmax_disp_lme2,tmax_disp_lme3,tmax_disp_lme4)

tmax_sens_lme4 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+year+year:seas+year:I(elev/1000)+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme4)
AIC(tmax_sens_lme1,tmax_sens_lme2,tmax_sens_lme3,tmax_sens_lme4)

##Adds a year:coldair interaction.
tavg_disp_lme5 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:cair_reg_9m+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme5)
AIC(tavg_disp_lme1,tavg_disp_lme2,tavg_disp_lme3,tavg_disp_lme4,tavg_disp_lme5)

tavg_sens_lme5 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:cair_reg_9m+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme5)
AIC(tavg_sens_lme1,tavg_sens_lme2,tavg_sens_lme3,tavg_sens_lme4,tavg_sens_lme5)

tmin_disp_lme5 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:cair_reg_9m+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme5)
AIC(tmin_disp_lme1,tmin_disp_lme2,tmin_disp_lme3,tmin_disp_lme4,tmin_disp_lme5)

tmin_sens_lme5 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)
                       +I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:cair_reg_9m+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme5)
AIC(tmin_sens_lme1,tmin_sens_lme2,tmin_sens_lme3,tmin_sens_lme4,tmin_sens_lme5)

tmax_disp_lme5 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:cair_reg_9m+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme5)
AIC(tmax_disp_lme1,tmax_disp_lme2,tmax_disp_lme3,tmax_disp_lme4,tmax_disp_lme5)

tmax_sens_lme5 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:cair_reg_9m+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme5)
AIC(tmax_sens_lme1,tmax_sens_lme2,tmax_sens_lme3,tmax_sens_lme4,tmax_sens_lme5)

##Substitutes a year:cancov interaction.
tavg_disp_lme6 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme6)
AIC(tavg_disp_lme1,tavg_disp_lme2,tavg_disp_lme3,tavg_disp_lme4,tavg_disp_lme5,tavg_disp_lme6)

tavg_sens_lme6 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme5)
AIC(tavg_sens_lme1,tavg_sens_lme2,tavg_sens_lme3,tavg_sens_lme4,tavg_sens_lme5,tavg_sens_lme6)

tmin_disp_lme6 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme6)
AIC(tmin_disp_lme1,tmin_disp_lme2,tmin_disp_lme3,tmin_disp_lme4,tmin_disp_lme5,tmin_disp_lme6)

tmin_sens_lme6 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)
                       +I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme6)
AIC(tmin_sens_lme1,tmin_sens_lme2,tmin_sens_lme3,tmin_sens_lme4,tmin_sens_lme5,tmin_sens_lme6)

tmax_disp_lme6 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme5)
AIC(tmax_disp_lme1,tmax_disp_lme2,tmax_disp_lme3,tmax_disp_lme4,tmax_disp_lme5,tmax_disp_lme6)

tmax_sens_lme6 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme5)
AIC(tmax_sens_lme1,tmax_sens_lme2,tmax_sens_lme3,tmax_sens_lme4,tmax_sens_lme5,tmax_sens_lme6)

##Substitutes a year:cancov interaction.
tavg_disp_lme6 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme6)
AIC(tavg_disp_lme1,tavg_disp_lme2,tavg_disp_lme3,tavg_disp_lme4,tavg_disp_lme5,tavg_disp_lme6)

tavg_sens_lme6 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme6)
AIC(tavg_sens_lme1,tavg_sens_lme2,tavg_sens_lme3,tavg_sens_lme4,tavg_sens_lme5,tavg_sens_lme6)

tmin_disp_lme6 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme6)
AIC(tmin_disp_lme1,tmin_disp_lme2,tmin_disp_lme3,tmin_disp_lme4,tmin_disp_lme5,tmin_disp_lme6)

tmin_sens_lme6 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)
                       +I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme6)
AIC(tmin_sens_lme1,tmin_sens_lme2,tmin_sens_lme3,tmin_sens_lme4,tmin_sens_lme5,tmin_sens_lme6)

tmax_disp_lme6 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme6)
AIC(tmax_disp_lme1,tmax_disp_lme2,tmax_disp_lme3,tmax_disp_lme4,tmax_disp_lme5,tmax_disp_lme6)

tmax_sens_lme6 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme6)
AIC(tmax_sens_lme1,tmax_sens_lme2,tmax_sens_lme3,tmax_sens_lme4,tmax_sens_lme5,tmax_sens_lme6)

##Substitutes a year:srad interaction.
tavg_disp_lme7 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(srad_seas/1e5)+year:I(srad_seas/1e5):seas+
                         cair_reg_9m:I(elev/1000)+(1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_disp_lme7)
AIC(tavg_disp_lme1,tavg_disp_lme2,tavg_disp_lme3,tavg_disp_lme4,tavg_disp_lme5,
    tavg_disp_lme6,tavg_disp_lme7)
best_mod_tavg_disp <- tavg_disp_lme7

tavg_sens_lme7 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+cair_reg_9m:I(elev/1000)+
                         year:I(X_UTM/1000)+
                         (1|site) + (1|site_name),
                       data=tavg_meta_sub)
summary(tavg_sens_lme7)
AIC(tavg_sens_lme1,tavg_sens_lme2,tavg_sens_lme3,tavg_sens_lme4,tavg_sens_lme5,
    tavg_sens_lme6,tavg_sens_lme7)
best_mod_tavg_sens <- tavg_sens_lme3

tmin_disp_lme7 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         year:I(X_UTM/1000):seas+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_disp_lme7)
AIC(tmin_disp_lme1,tmin_disp_lme2,tmin_disp_lme3,tmin_disp_lme4,tmin_disp_lme5,
    tmin_disp_lme6,tmin_disp_lme7)
best_mod_tmin_disp <- tmin_disp_lme7
  
tmin_sens_lme7 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)
                       +I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         year:I(X_UTM/1000)+year:I(X_UTM/1000):seas+
                         (1|site) + (1|site_name),
                       data=tmin_meta_sub)
summary(tmin_sens_lme7)
AIC(tmin_sens_lme1,tmin_sens_lme2,tmin_sens_lme3,tmin_sens_lme4,tmin_sens_lme5,
    tmin_sens_lme6,tmin_sens_lme7)
best_mod_tmin_sens <- tmin_sens_lme4

tmax_disp_lme7 <- lmer(disp~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+I(srad_seas/1e5):I(X_UTM/1000)+I(srad_seas/1e5):I(Y_UTM/1000)+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         year:I(reg_cancov_243m/100):seas+
                         (1|site) + (1|site_name),
                       data=tmax_meta_sub)
summary(tmax_disp_lme7)
AIC(tmax_disp_lme1,tmax_disp_lme2,tmax_disp_lme3,tmax_disp_lme4,tmax_disp_lme5,
    tmax_disp_lme6,tmax_disp_lme7)
best_mod_tmax_disp <- tmax_disp_lme2

tmax_sens_lme7 <- lmer(sens~seas+I(elev/1000)+cair_reg_9m+I(reg_cancov_243m/100)+
                         I(srad_seas/1e5)+I(srad_seas/1e5):seas+
                         I(X_UTM/1000)+I(Y_UTM/1000)+cair_reg_9m:seas+I(elev/1000):seas+type+
                         year+year:seas+year:I(elev/1000)+year:I(reg_cancov_243m/100)+
                         year:I(reg_cancov_243m/100):seas+
                         (1|site)+(1|site_name),
                       data=tmax_meta_sub)
summary(tmax_sens_lme7)
AIC(tmax_sens_lme1,tmax_sens_lme2,tmax_sens_lme3,tmax_sens_lme4,tmax_sens_lme5,
    tmax_sens_lme6,tmax_sens_lme7)
best_mod_tmax_sens <- tmax_sens_lme1

formula(best_mod_tavg_disp)
formula(best_mod_tavg_sens)
formula(best_mod_tmax_disp)
formula(best_mod_tmax_sens)
formula(best_mod_tmin_disp)
formula(best_mod_tmin_sens)


##Picks best model by AIC and predicts from it to make graphs.
elevs <- c(seq(500,2200,by=50))
srads_DJF <- quantile(tmax_meta_sub$srad_seas_DJF,probs=c(0.05,0.95))
srads_MAM <- quantile(tmax_meta_sub$srad_seas_MAM,probs=c(0.05,0.95))
srads_JJA <- quantile(tmax_meta_sub$srad_seas_JJA,probs=c(0.05,0.95))
srads_SON <- quantile(tmax_meta_sub$srad_seas_SON,probs=c(0.05,0.95))
canopies <- c(0.05,0.5,0.95)
cairs <- quantile(tmax_meta_sub$cair_reg_9m,probs=c(0.05,0.25,0.5,0.75,0.95))
xs <- quantile(tmax_meta_sub$X_UTM)[c(2,3,4)]
ys <- quantile(tmax_meta_sub$Y_UTM)[c(2,3,4)]
seas <- c("DJF","MAM","JJA","SON")
years <- unique(tmax_meta_sub$year)

pred_data_MAM_buff <- expand.grid(elev=elevs,
                                  srad_seas=srads_MAM[1],
                                  seas="MAM",
                                  reg_cancov_243m=canopies[3],
                                  X_UTM=xs[2],
                                  Y_UTM=ys[2],
                                  cair_reg_9m=cairs[5],
                                  year=years,
                                  type="DataLogger")
pred_data_MAM_buff$SiteType <- "Shaded Deep-forest Cove"
pred_data_MAM_sens <- expand.grid(elev=elevs,
                                  srad_seas=srads_MAM[2],
                                  seas="MAM",
                                  reg_cancov_243m=canopies[1],
                                  X_UTM=xs[2],
                                  Y_UTM=ys[2],
                                  cair_reg_9m=cairs[1],
                                  year=years,
                                  type="DataLogger")
pred_data_MAM_sens$SiteType <- "Sunny Open-canopy Ridge"
pred_data_MAM <- rbind(pred_data_MAM_buff,pred_data_MAM_sens)

boot_pred_MAM_fun <- function(x){
  predict(x,newdata=pred_data_MAM,re.form=NA)}
quantfun <- function(x) {quantile(x,probs=c(0.025,0.975))}

pred_data_MAM$tavg_pred <- predict(best_mod_tavg_disp,newdata=pred_data_MAM,re.form=NA)
pred_ci_MAM_disp <- t(apply(bootMer(best_mod_tavg_disp,nsim=500,FUN=boot_pred_MAM_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_MAM$tavg_lwr <- pred_ci_MAM_disp[,1]
pred_data_MAM$tavg_upr <- pred_ci_MAM_disp[,2]
pred_data_MAM$tmin_pred <- predict(best_mod_tmin_disp,newdata=pred_data_MAM,re.form=NA)
pred_ci_MAM_disp <- t(apply(bootMer(best_mod_tmin_disp,nsim=500,FUN=boot_pred_MAM_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_MAM$tmin_lwr <- pred_ci_MAM_disp[,1]
pred_data_MAM$tmin_upr <- pred_ci_MAM_disp[,2]
pred_data_MAM$tmax_pred <- predict(best_mod_tmax_disp,newdata=pred_data_MAM,re.form=NA)
pred_ci_MAM_disp <- t(apply(bootMer(best_mod_tmax_disp,nsim=500,FUN=boot_pred_MAM_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_MAM$tmax_lwr <- pred_ci_MAM_disp[,1]
pred_data_MAM$tmax_upr <- pred_ci_MAM_disp[,2]


pred_data_JJA_buff <- expand.grid(elev=elevs,
                             srad_seas=srads_JJA[1],
                             seas="JJA",
                             reg_cancov_243m=canopies[3],
                             X_UTM=xs[2],
                             Y_UTM=ys[2],
                             cair_reg_9m=cairs[5],
                             year=years,
                             type="DataLogger")
pred_data_JJA_buff$SiteType <- "Shaded Deep-forest Cove"
pred_data_JJA_sens <- expand.grid(elev=elevs,
                                  srad_seas=srads_JJA[2],
                                  seas="JJA",
                                  reg_cancov_243m=canopies[1],
                                  X_UTM=xs[2],
                                  Y_UTM=ys[2],
                                  cair_reg_9m=cairs[1],
                                  year=years,
                                  type="DataLogger")
pred_data_JJA_sens$SiteType <- "Sunny Open-canopy Ridge"
pred_data_JJA <- rbind(pred_data_JJA_buff,pred_data_JJA_sens)

boot_pred_JJA_fun <- function(x){
  predict(x,newdata=pred_data_JJA,re.form=NA)}
quantfun <- function(x) {quantile(x,probs=c(0.025,0.975))}

pred_data_JJA$tavg_pred <- predict(best_mod_tavg_disp,newdata=pred_data_JJA,re.form=NA)
pred_ci_JJA_disp <- t(apply(bootMer(best_mod_tavg_disp,nsim=500,FUN=boot_pred_JJA_fun,
                                type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_JJA$tavg_lwr <- pred_ci_JJA_disp[,1]
pred_data_JJA$tavg_upr <- pred_ci_JJA_disp[,2]
pred_data_JJA$tmin_pred <- predict(best_mod_tmin_disp,newdata=pred_data_JJA,re.form=NA)
pred_ci_JJA_disp <- t(apply(bootMer(best_mod_tmin_disp,nsim=500,FUN=boot_pred_JJA_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_JJA$tmin_lwr <- pred_ci_JJA_disp[,1]
pred_data_JJA$tmin_upr <- pred_ci_JJA_disp[,2]
pred_data_JJA$tmax_pred <- predict(best_mod_tmax_disp,newdata=pred_data_JJA,re.form=NA)
pred_ci_JJA_disp <- t(apply(bootMer(best_mod_tmax_disp,nsim=500,FUN=boot_pred_JJA_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_JJA$tmax_lwr <- pred_ci_JJA_disp[,1]
pred_data_JJA$tmax_upr <- pred_ci_JJA_disp[,2]

pred_data_SON_buff <- expand.grid(elev=elevs,
                                  srad_seas=srads_SON[1],
                                  seas="SON",
                                  reg_cancov_243m=canopies[3],
                                  X_UTM=xs[2],
                                  Y_UTM=ys[2],
                                  cair_reg_9m=cairs[5],
                                  year=years,
                                  type="DataLogger")
pred_data_SON_buff$SiteType <- "Shaded Deep-forest Cove"
pred_data_SON_sens <- expand.grid(elev=elevs,
                                  srad_seas=srads_SON[2],
                                  seas="SON",
                                  reg_cancov_243m=canopies[1],
                                  X_UTM=xs[2],
                                  Y_UTM=ys[2],
                                  cair_reg_9m=cairs[1],
                                  year=years,
                                  type="DataLogger")
pred_data_SON_sens$SiteType <- "Sunny Open-canopy Ridge"
pred_data_SON <- rbind(pred_data_SON_buff,pred_data_SON_sens)

boot_pred_SON_fun <- function(x){
  predict(x,newdata=pred_data_SON,re.form=NA)}
quantfun <- function(x) {quantile(x,probs=c(0.025,0.975))}

pred_data_SON$tavg_pred <- predict(best_mod_tavg_disp,newdata=pred_data_SON,re.form=NA)
pred_ci_SON_disp <- t(apply(bootMer(best_mod_tavg_disp,nsim=500,FUN=boot_pred_SON_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_SON$tavg_lwr <- pred_ci_SON_disp[,1]
pred_data_SON$tavg_upr <- pred_ci_SON_disp[,2]
pred_data_SON$tmin_pred <- predict(best_mod_tmin_disp,newdata=pred_data_SON,re.form=NA)
pred_ci_SON_disp <- t(apply(bootMer(best_mod_tmin_disp,nsim=500,FUN=boot_pred_SON_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_SON$tmin_lwr <- pred_ci_SON_disp[,1]
pred_data_SON$tmin_upr <- pred_ci_SON_disp[,2]
pred_data_SON$tmax_pred <- predict(best_mod_tmax_disp,newdata=pred_data_SON,re.form=NA)
pred_ci_SON_disp <- t(apply(bootMer(best_mod_tmax_disp,nsim=500,FUN=boot_pred_SON_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_SON$tmax_lwr <- pred_ci_SON_disp[,1]
pred_data_SON$tmax_upr <- pred_ci_SON_disp[,2]

pred_data_DJF_buff <- expand.grid(elev=elevs,
                                  srad_seas=srads_DJF[1],
                                  seas="DJF",
                                  reg_cancov_243m=canopies[3],
                                  X_UTM=xs[2],
                                  Y_UTM=ys[2],
                                  cair_reg_9m=cairs[5],
                                  year=years,
                                  type="DataLogger")
pred_data_DJF_buff$SiteType <- "Shaded Deep-forest Cove"
pred_data_DJF_sens <- expand.grid(elev=elevs,
                                  srad_seas=srads_DJF[2],
                                  seas="DJF",
                                  reg_cancov_243m=canopies[1],
                                  X_UTM=xs[2],
                                  Y_UTM=ys[2],
                                  cair_reg_9m=cairs[1],
                                  year=years,
                                  type="DataLogger")
pred_data_DJF_sens$SiteType <- "Sunny Open-canopy Ridge"
pred_data_DJF <- rbind(pred_data_DJF_buff,pred_data_DJF_sens)

boot_pred_DJF_fun <- function(x){
  predict(x,newdata=pred_data_DJF,re.form=NA)}
quantfun <- function(x) {quantile(x,probs=c(0.025,0.975))}

pred_data_DJF$tavg_pred <- predict(best_mod_tavg_disp,newdata=pred_data_DJF,re.form=NA)
pred_ci_DJF_disp <- t(apply(bootMer(best_mod_tavg_disp,nsim=500,FUN=boot_pred_DJF_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_DJF$tavg_lwr <- pred_ci_DJF_disp[,1]
pred_data_DJF$tavg_upr <- pred_ci_DJF_disp[,2]
pred_data_DJF$tmin_pred <- predict(best_mod_tmin_disp,newdata=pred_data_DJF,re.form=NA)
pred_ci_DJF_disp <- t(apply(bootMer(best_mod_tmin_disp,nsim=500,FUN=boot_pred_DJF_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_DJF$tmin_lwr <- pred_ci_DJF_disp[,1]
pred_data_DJF$tmin_upr <- pred_ci_DJF_disp[,2]
pred_data_DJF$tmax_pred <- predict(best_mod_tmax_disp,newdata=pred_data_DJF,re.form=NA)
pred_ci_DJF_disp <- t(apply(bootMer(best_mod_tmax_disp,nsim=500,FUN=boot_pred_DJF_fun,
                                    type="parametric",re.form=NA)$t,
                            FUN=quantfun,MARGIN=2))
pred_data_DJF$tmax_lwr <- pred_ci_DJF_disp[,1]
pred_data_DJF$tmax_upr <- pred_ci_DJF_disp[,2]

pred_data_allseas <- rbind(pred_data_DJF,
                           pred_data_MAM,
                           pred_data_JJA,
                           pred_data_SON)
##Plots predictions.
p1 <- ggplot(pred_data_allseas)+
  geom_hline(aes(yintercept=0),alpha=0.5)+
  geom_line(aes(x=elev,y=tmin_pred,color=SiteType,linetype=SiteType,fill=year))+
  geom_ribbon(aes(x=elev,ymin=tmin_lwr,ymax=tmin_upr,
                  fill=SiteType,alpha=year))+
  scale_alpha_manual(values=c(0.1,0.1,0.1,0.1,0.1,0.1))+
  scale_linetype_manual(values=c(1,1))+
  scale_fill_grey(start=0.2,end=0.6)+
  scale_color_grey(start=0.2,end=0.6)+
  ylim(c(-11,5))+
  ylab("Tmin")+
  xlim(c(500,2200))+
  xlab("")+
  facet_grid(facets=.~seas)+
  theme_bw()+
  theme(legend.position="none",legend.direction="horizontal",
        panel.grid=element_blank())
##Plots predictions.
p2 <- ggplot(pred_data_allseas)+
  geom_hline(aes(yintercept=0),alpha=0.5)+
  geom_line(aes(x=elev,y=tmax_pred,color=SiteType,linetype=SiteType,fill=year))+
  geom_ribbon(aes(x=elev,ymin=tmax_lwr,ymax=tmax_upr,
                  fill=SiteType,alpha=year))+
  scale_alpha_manual(values=c(0.1,0.1,0.1,0.1,0.1,0.1))+
  scale_linetype_manual(values=c(1,1))+
  scale_fill_grey(start=0.2,end=0.6)+
  scale_color_grey(start=0.2,end=0.6)+
  ylab("Tmax")+
  ylim(c(-11,5))+
  xlab("")+
  xlim(c(500,2200))+
  facet_grid(facets=.~seas)+
  theme_bw()+
  theme(legend.position="none",legend.direction="horizontal",
        panel.grid=element_blank())
p3 <- ggplot(pred_data_allseas)+
  geom_hline(aes(yintercept=0),alpha=0.5)+
  geom_line(aes(x=elev,y=tavg_pred,color=SiteType,linetype=SiteType,fill=year))+
  geom_ribbon(aes(x=elev,ymin=tavg_lwr,ymax=tavg_upr,
                  fill=SiteType,alpha=year))+
  scale_alpha_manual(values=c(0.1,0.1,0.1,0.1,0.1,0.1),guide=FALSE)+
  scale_linetype_manual(values=c(1,1),guide=FALSE)+
  scale_fill_grey(start=0.2,end=0.6,guide=FALSE)+
  scale_color_grey(start=0.2,end=0.6)+
  ylab("Tavg")+
  ylim(c(-11,5))+
  xlim(c(500,2200))+
  xlab("Elevation (m)")+
  facet_grid(facets=.~seas)+
  theme_bw()+
  theme(legend.position="bottom",legend.direction="horizontal",
        panel.grid=element_blank())

require(gridExtra)
pdf("./results/reg_microclim_yearly_variation.pdf",width=7,height=6)
grid.arrange(p1,p2,p3,nrow=3,heights=c(1,1,1.5))
dev.off()


####Maps the residuals of the best model.

##puts the best models in a list
best_mods <- list(best_mod_tavg_disp,best_mod_tavg_sens,
                  best_mod_tmin_disp,best_mod_tavg_sens,
                  best_mod_tmax_disp,best_mod_tmax_sens)

resid_maps <- data.frame(X=tmax_meta_sub$X_UTM,
                         Y=tmax_meta_sub$Y_UTM,
                         resid=residuals(best_mods[[5]]))
ggplot(resid_maps)+
  geom_point(aes(x=X,y=Y,color=resid))+
  scale_color_gradient2()+
  theme_bw()

####Predicts regional microclimate grids based on the best model.####
library(raster)
elev_9m <- raster("/Volumes/ib_working/GIS/elev_NED_9m.tif")
cold_9m <- raster("/Volumes/ib_working/GIS/coldair_index_regional_9m.tif")
reg_cancov_27m <- raster("/Volumes/ib_working/GIS/reg_canopy_focal27m.tif")
reg_cancov_81m <- raster("/Volumes/ib_working/GIS/reg_canopy_focal81m.tif")
reg_cancov_243m <- raster("/Volumes/ib_working/GIS/reg_canopy_focal243m.tif")
srad_reg_DJF <- raster("/Volumes/ib_working/GIS/srad_seas_DJF.tiff")
srad_reg_MAM <- raster("/Volumes/ib_working/GIS/srad_seas_MAM.tiff")
srad_reg_JJA <- raster("/Volumes/ib_working/GIS/srad_seas_JJA.tiff")
srad_reg_SON <- raster("/Volumes/ib_working/GIS/srad_seas_SON.tiff")
X_UTM <- raster("/Volumes/ib_working/GIS/reg_UTM_X.tif")
Y_UTM <- raster("/Volumes/ib_working/GIS/reg_UTM_Y.tif")
reg_stack <- stack(elev_9m,cold_9m,reg_cancov_243m,
                   srad_reg_DJF,srad_reg_MAM,srad_reg_JJA,srad_reg_SON,
                   X_UTM,Y_UTM)
reg_brick <- brick(reg_stack,filename=c("~/GIS/reg_mclim_preds.grd"),datatype="FLT4S")
names(reg_brick) <- c("elev","cair_reg_9m","reg_cancov_243m","srad_seas","srad_seas_MAM",
                     "srad_seas_JJA","srad_seas_SON","X_UTM","Y_UTM")



test_ext <- extent(c(reg_brick@extent@xmin+16000,
                     reg_brick@extent@xmin+27000,
                     reg_brick@extent@ymin+7000,
                     reg_brick@extent@ymin+21000))
test_env <- crop(reg_brick,test_ext,filename="~/GIS/reg_microclim_preds_9m_test.grd",overwrite=TRUE)
names(test_env) <- c("elev","cair_reg_9m","reg_cancov_243m","srad_seas","srad_seas_MAM",
                      "srad_seas_JJA","srad_seas_SON","X_UTM","Y_UTM")

##Tests the raster predictions.

test_pred1 <- raster::predict(test_env,best_mods[[1]],re.form=~0,
                              const=data.frame(seas="DJF",
                                               year="2009-10",
                                               type="DataLogger"))

test_pred2 <- raster::predict(test_env,best_mods[[1]],re.form=~0,
                             const=data.frame(seas="DJF",
                                              year="2011-12",
                                              type="DataLogger"))
test_pred3 <- raster::predict(test_env,best_mods[[1]],re.form=~0,
                              const=data.frame(seas="DJF",
                                               year="2009-10",
                                               type="DataLogger"))
test_pred4 <- raster::predict(test_env,best_mods[[1]],re.form=~0,
                              const=data.frame(seas="DJF",
                                               year="2014-15",
                                               type="DataLogger"))
plot(stack(test_pred1,test_pred2,test_pred3,test_pred4),zlim=c(0,5))

##Tests the full-extent predictions.
test_full <- raster::predict(reg_brick,best_mods[[1]],re.form=~0,
                             const=data.frame(seas="DJF",
                                              year="2014-15",
                                              type="DataLogger"),
                             filename=c("/Volumes/ib_working/GIS/test_DJF_full.tif"))
