##Script to build statistical models of snow duration.
##Author:Ian Breckheimer

####Sets up the workspace####

## Loads required packages
library(lme4)
library(ggmap)
library(raster)
library(ggplot2)
library(gstat)
library(automap)
library(dplyr)

##Logit and inverse logit function.
logit <- function(x){log(x/(1-x))}
inv.logit <- function(x){exp(x)/(1+exp(x))}

## Sets working directory
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/cleaned")

##Loads snow duration data.
snow_output <- read.csv("microclimate_snow_coords_final_10_19_15.csv")
snow_output$water_year <- as.factor(snow_output$water_year)

####Exploratory Analysis####
##Quick plots to check data.
qplot(snow_cover_days,elevation,data=snow_output,facets=water_year~.,color=relev_30m,xlim=c(0,365))+
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y~poly(x,2))+
  labs(x="Snow Duration (Days)",
       y="Elevation (m)",
       color="Topo\nExposure",
       title="Water Year")+
  theme_bw()

##Quick plots to check data.
snow_output$yearfact <- as.factor(snow_output$water_year)
pdf("~/Dropbox/MORA_phenology/figs/snow_2009_2015.pdf",width=6,height=4)
ggplot(data=snow_output)+
  geom_point(aes(y=snow_diss_DOY,x=elevation,color=yearfact),alpha=0)+
  geom_smooth(aes(y=snow_diss_DOY,x=elevation,color=yearfact,linetype=yearfact),
              method = "lm", se=TRUE, formula = y~poly(x,2))+
  scale_color_brewer("Year",palette = "Dark2")+
  scale_linetype_discrete("Year")+
  labs(y="Snow Dissapearance Date (Julian Day)",
       x="Elevation (m)")+
  xlim(c(600,2200))+
  ylim(c(0,250))+
  theme_bw()+
  theme(panel.grid=element_blank())
dev.off()



####Preliminary statistical models of SDD####

##Only 2 obs for 2007 so remove.
snow_output_dat <- subset(snow_output,water_year %in% c("2009","2010","2011","2012","2013","2014","2015"))

####Replot with 2007 removed.
##Quick plots to check data.

##Color scale consistent with map.
snow_cols <- c(rgb(250/255,8/255,8/255),
               rgb(241/255,165/255,191/255),
               rgb(123/255,172/255,30/255),
               rgb(50/255,98/255,132/255),
               rgb(136/255,7/255,129/255),
               rgb(0,0,0))

# ##Melt dates for paradise SNOTEL 2009 - 2015
para_dates <- c(202,205,240,209,199,204,151)
para_elevs <- rep(1563,7)
#
##Appends those values to data frame.
snow_output_para <- snow_output_dat[1:7,]
snow_output_para[,] <- rep(NA,203)
snow_output_para$snow_diss_DOY <- para_dates
snow_output_para$elevation <- para_elevs
snow_output_para$study <- rep("Para_SNOTEL",7)
snow_output_para$site_name <- rep("Para_SNOTEL",7)
snow_output_para$water_year <- 2009:2015
snow_output_para$PARA <- TRUE

snow_output_dat$PARA <- FALSE
snow_output_dat_para <- rbind(snow_output_dat,snow_output_para)


cairo_pdf("sdd_plot_allyears_2015.pdf",width=2.5,height=5,family="CMU Serif")
ggplot(data=snow_output_dat_para)+
  facet_grid(water_year~.)+
  geom_point(aes(x=snow_diss_DOY,y=elevation,color=snow_out,size=PARA,shape=PARA))+
  geom_smooth(aes(x=snow_diss_DOY,y=elevation),method = "lm", se=TRUE, color="black", formula = y~poly(x,2))+
  scale_colour_gradientn(colours=snow_cols, values = NULL,
                         space = "Lab", na.value = "grey50",
                         guide = "colourbar")+
  scale_size_manual(values=c(1,4))+
  scale_x_continuous(breaks=c(32,91,152,213), 
                     labels=c("Feb.","April", "June","August"))+
  labs(x="Snow Melt Date",
       y="Elevation (m)",
       color="Canopy\nCover (%)",
       title="")+
  theme_bw()+
  theme(legend.position="none",
        panel.border=element_rect(fill=NULL,color="black",size=1),
        strip.background=element_rect(fill="white",color="white"))
dev.off()

snowmod <- lm(snow_diss_DOY~elevation+water_year,data=snow_output_dat)
summary(snowmod)
snowmod2 <- lm(snow_diss_DOY~elevation*water_year,data=snow_output_dat)
summary(snowmod2)
snowmod3 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2),data=snow_output_dat)
summary(snowmod3)
snowmod4 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM + Y_UTM,data=snow_output_dat)
summary(snowmod4)
snowmod5 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM*Y_UTM,data=snow_output_dat)
summary(snowmod5)
snowmod6 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM*Y_UTM+canopy_pct,data=snow_output_dat)
summary(snowmod6)
snowmod7 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope,data=snow_output_dat)
summary(snowmod7)
snowmod8 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m,data=snow_output_dat)
summary(snowmod8)
snowmod9 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m*elevation+srad,data=snow_output_dat)
summary(snowmod9)
snowmod10 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m*elevation+srad+stream_dist,data=snow_output_dat)
summary(snowmod10)
snowmod11 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)*water_year+X_UTM*Y_UTM+X_UTM*water_year+canopy_pct+slope+relev_30m*elevation+srad+stream_dist,data=snow_output_dat)
summary(snowmod11)
AIC(snowmod,snowmod2,snowmod3,snowmod4,snowmod5,snowmod6,snowmod7,snowmod8,snowmod9,snowmod10,snowmod11)

####Examines spatial autocorrelation in the residuals####

## Subsets the data
snow_spat <- subset(snow_output,water_year %in% c("2009","2010","2011","2012","2013","2014","2015"))

## Residuals from best model
snow_spat$resid <- snowmod10$residuals

## Converts coordinates to km to avoid numerical problems with fitting variograms
obs_coords <- as.matrix(snow_spat[,c("X_UTM","Y_UTM")])/1000

## Converts data to a spatial points data frame.
coordinates(snow_spat) <- obs_coords

## Computes the empirical variogram on the raw data and residuals, 
## then fits a variogram model to the data.
snow_vg <- variogram(snow_diss_DOY~1,data=snow_spat,
                  width=0.1,cutoff=10)
snow_fit <- autofitVariogram(snow_diss_DOY~1,input_data=snow_spat,model="Exp",
                             fix.values=c(0,NA,NA),width=0.05,cutoff=10)

resid_vg <- variogram(resid~1,data=snow_spat,
                      width=0.1,cutoff=10)
resid_fit <- autofitVariogram(resid~1,input_data=snow_spat,model="Exp",
              fix.values=c(0,NA,NA),width=0.05,cutoff=10)
## Plots the variogram models.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/maps")
svg("snow_variogram.svg",width=4,height=4,family='sans')
par(mfrow=c(2,1),mar=c(0,3,0,0),oma=c(4,1,1,1))
plot(snow_vg,ylim=c(0,1600),pch=19,col=rgb(0,0.7,0.9,0.5),cex=0.5,xaxt='n')
lines.variomodel(snow_fit,lwd=1.5)
text(2,50,labels="Raw Data",pos=2)
plot(resid_vg,pch=19,col=rgb(0,0.7,0.2,0.5),cex=0.5,ylim=c(0,1600))
lines.variomodel(resid_fit,lwd=1.5)
text(2,400,labels="Model Residuals",pos=2)
mtext("Distance (km)",side=1,outer=T,padj = 4)
mtext("Semivariance",side=2,outer=T,padj = 1)
dev.off()

####Mixed-effects model to estimate snow melt lapse rate.
snowmod_avg <- lmer(snow_diss_DOY~elevation+I(scale(X_UTM,scale=FALSE))+I(scale(Y_UTM,scale=FALSE))+water_year + (1|site_name),
                    data=snow_output_dat,REML=FALSE)
summary(snowmod_avg)

####Mixed-effects models to deal with non-independence by location####
snowmod_lme <- lmer(snow_diss_DOY~elevation + water_year + (1 |site_name),
                    data=snow_output_dat,REML=FALSE)
summary(snowmod_lme)
snowmod_lme2 <- lmer(snow_diss_DOY~elevation + (1 |site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme2)
snowmod_lme3 <- lmer(snow_diss_DOY~elevation + canopy_pct + (1 |site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme3)
snowmod_lme4 <- lmer(snow_diss_DOY~elevation + canopy_pct + I(Y_UTM/1000) + (1 |site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme4)
snowmod_lme5 <- lmer(snow_diss_DOY~elevation + canopy_pct + relev_30m + I(Y_UTM/1000) + (1 |site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme5)
snowmod_lme6 <- lmer(snow_diss_DOY~water_year+elevation + I(elevation^2) + canopy_pct + slope + relev_30m + I(Y_UTM/1000) + I(X_UTM/1000) + (1 |site_name),
                     data=snow_output_dat,REML=FALSE,control=lmerControl(optCtrl=list(maxfun=20000)))
summary(snowmod_lme6)
snowmod_lme7 <- lmer(snow_diss_DOY~elevation*water_year+I(elevation^2)*water_year+I(X_UTM/1000)*water_year + canopy_pct + slope +
                       relev_30m*elevation+I(srad_canopy_30m/1000) + (1 |site_name),
                     data=snow_output_dat,REML=FALSE,control=lmerControl(optCtrl=list(maxfun=20000)))
summary(snowmod_lme7)
resid <- residuals(snowmod_lme7)
plot(density(resid))
AIC(snowmod_lme,snowmod_lme2,snowmod_lme3,snowmod_lme4,snowmod_lme5,snowmod_lme6,snowmod_lme7)

##Best lme model with REML and minus outliers for estimating parameters
snow_output_drop <- snow_output_dat[-c(which(residuals(snowmod_lme7) > 40| residuals(snowmod_lme7) < -40)),]
snowmod_lme8 <- lmer(snow_diss_DOY~elevation*water_year+I(elevation^2)*water_year+I(X_UTM/1000)*water_year + canopy_pct+
                       slope + relev_30m*elevation + I(srad/1000) + (1|site_name),
                     data=snow_output_drop,REML=TRUE,control=lmerControl(optCtrl=list(maxfun=40000)))
summary(snowmod_lme8)

##Computes bootstrap confidence intervals for the parameters.
lme8_coef <- fixef(snowmod_lme8)
lme8_conf <- confint(snowmod_lme8,method="boot",nsim=1000)
lme8_table <- cbind(lme8_coef,data.frame(lme8_conf)[-c(1:2),])
write.csv(lme8_table,"/Users/ian/code/MORA_pheno/scratch/snow_lme_params.csv",row.names=TRUE)

##Plots residuals of the final model.
resid8 <- residuals(snowmod_lme8)
predict8 <- predict(snowmod_lme8,re.form=~0)
snow_output_drop$lme_pred <- predict8
coordinates(snow_output_drop) <- ~UTMX+UTMY

pdf("/Users/ian/code/MORA_pheno/figs/Snow_model_validation.pdf",width=8,height=5)
par(mfrow=c(1,2),mar=c(4,4,2,0))
plot(density(resid8),main="",xlab="SDD Residual (Days)",
     ylab="Density")
plot(lme_pred~snow_diss_DOY,data=snow_output_drop,pch=20,cex=0.3,col=factor(snow_output_drop$water_year),
     xlab="Observed SDD (Day of Year)",
     ylab="Predicted SDD (Day of Year)",
     xlim=c(0,250),
     ylim=c(0,250))
abline(0,1,lty=2)
dev.off()

##Puts the lme residuals on a map
snow_output_drop$lmeresid <- residuals(snowmod_lme8)
map <- get_map(location=c(-121.75,46.85),zoom=11,maptype="satellite",source="google")
ggmap(map)+
  geom_point(aes(x = jitter(longitude,amount=0.005), 
                 y = jitter(latitude,amount=0.005), 
                 color = lmeresid), 
             data = subset(snow_output_drop,lmeresid > 10 | lmeresid < -10) , alpha = .8)+
  scale_colour_gradientn(colours = c("#EDAAAA","#FFFFFF","#0800FF"),values=c(0,0.5,1))

##Plots a semivariogram of residuals.
## Converts coordinates to km to avoid numerical problems with fitting variograms
snow_spat <- snow_output_drop
obs_coords <- as.matrix(snow_output_drop[,c("X_UTM","Y_UTM")])/1000

## Converts data to a spatial points data frame.
## Converts coordinates to km to avoid numerical problems with fitting variograms
obs_coords <- as.matrix(snow_output_drop[,c("X_UTM","Y_UTM")])/1000

## Converts data to a spatial points data frame.
coordinates(snow_output_drop) <- obs_coords


## Computes the empirical variogram on the raw data and residuals, 
## then fits a variogram model to the data.

snow_vg_lme <- variog(coords=snow_output_drop@coords,data=snow_output_drop@data$snow_diss_DOY,option="bin",
                  breaks=seq(0,8,by=0.01),max.dist=2)
snow_fit_lme <- variofit(snow_vg_lme,ini.cov.pars=c(600,0.5),nugget=0,cov.model="exp",fix.nugget=FALSE)

resid_vg_lme <- variog(coords=snow_output_drop@coords,data=snow_output_drop@data$lmeresid,option="bin",
                   breaks=seq(0,8,by=0.01),max.dist=2)
resid_fit_lme <- variofit(resid_vg_lme,ini.cov.pars=c(200,0.6),nugget=0,cov.model="exp",fix.nugget=FALSE)

pdf("/Users/ian/code/MORA_pheno/figs/Snow_model_variogram.pdf",width=4,height=5,family='sans')
par(mfrow=c(2,1),mar=c(0,3,0,0),oma=c(4,1,2,2))
plot(snow_vg_lme,ylim=c(0,940),pch=19,col=rgb(0,0.7,0.9,0.5),cex=0.5,xaxt='t')
lines.variomodel(snow_fit_lme,lwd=1.5)
text(2,50,labels="Raw Data",pos=2)
plot(resid_vg_lme,pch=19,col=rgb(0,0.7,0.2,0.5),cex=0.5,ylim=c(0,940))
lines.variomodel(resid_fit_lme,lwd=1.5)
text(2,400,labels="Model Residuals",pos=2)
mtext("Distance (km)",side=1,outer=T,padj = 4)
mtext("Semivariance",side=2,outer=T,padj = 1)
dev.off()

####Generates snow predictions with uncertainty for all of the flickr points.####
flickr <- read.csv("~/code/MORA_pheno/data/MORA_flickr_metadata_all_2009_2015_covar_final.csv")
flickr$water_year <- as.factor(flickr$year)

flickr$sdd_pred <- predict(snowmod_lme8,flickr,re.form=~0)

##Bootstraps prediction intervals.
predFun <- function(fit) {
  predict(fit,flickr,re.form=~0)
}
fb <- bootMer(snowmod_lme8,nsim=1000,FUN=predFun,seed=101)

##Extracts quantiles from the bootstrap samples.
quantile_fun_lwr <- function(x) {stats::quantile(x,probs=c(0.025),na.rm=TRUE)}
quantile_fun_upr <- function(x) {stats::quantile(x,probs=c(0.975),na.rm=TRUE)}

flickr$sdd_lwr <- apply(fb$t,FUN = quantile_fun_lwr,MARGIN=2)
flickr$sdd_upr <- apply(fb$t,FUN = quantile_fun_upr,MARGIN=2)
flickr$sdd_range <- flickr$sdd_upr - flickr$sdd_lwr

write.csv(flickr,"~/code/MORA_pheno/data/MORA_flickr_metadata_all_2009_2015_covar_sdd_final.csv",
          row.names=FALSE)

####Tries fitting a GAMM.
library(gamm4)
gamm_mod1 <- gamm4(snow_diss_DOY~s(elevation,water_year,bs="fs",k=8)+I(X_UTM/1000)*water_year +I(Y_UTM/1000)*water_year + canopy_pct + slope  + 
                  relev_30m + relev_30m:elevation + I(srad/1000),family=gaussian(),
                 data=snow_output_dat, random = ~(1|site_name))
summary(gamm_mod1$mer)
summary(gamm_mod1$gam)
plot(gamm_mod1$gam)
####Constructs raster predictions for each year based on the best mixed-effects model####

# gf <- focalWeight(env3$elev, 30, "circle")
# env3$relev_30m <- env3$elev - focal(env3$elev,gf)
# env3$X_UTM <- raster('/Volumes/ib_working/GIS/MORA_UTMX_3m.tif')
# env3$Y_UTM <- raster('/Volumes/ib_working/GIS/MORA_UTMY_3m.tif')
# writeRaster(env3,"~/GIS/MORA_snow_env_pred2.grd",format="raster",overwrite=T)
# env_pred <- brick(env3[[c(3,4,6,7,14,15,17:21)]],filename="~/GIS/MORA_snow_pred.grd",
#                    format="raster",overwrite=T)
env_pred <- brick("/Volumes/ib_working/GIS/MORA_snow_pred.grd")

names(env_pred) <- c("canopy_pct","elevation","slope","srad","stream_dist","asp",
                     "cold_air_ind","srad_noc","relev_30m","X_UTM","Y_UTM")

## Subsets the data
snow_output_dat <- subset(snow_output,water_year %in% c("2009","2010","2011",
                                                        "2012","2013","2014","2015"))
par(mfrow=c(1,1),mar=c(4,4,1,1))
pred <- predict(snowmod_lme8,re.form=NA)
plot(snow_output_drop$snow_diss_DOY~pred,xlab="Observed",ylab="Predicted")
abline(0,1,lty=2)

##Computes the average absolute error.
mean(abs(snow_output_drop$snow_diss_DOY - pred))

####Computes raster predictions for 2009-2015 based on the best mixed-effects model.

##Tweaks options for speed
rasterOptions(maxmemory=5e+08,chunksize=1e+07,timer=TRUE)

##Tests the model on a small extent.
test_ext <- extent(c(env_pred@extent@xmin+15000,
                     env_pred@extent@xmin+18000,
                     env_pred@extent@ymin+7000,
                     env_pred@extent@ymin+11000))
test_env <- crop(env_pred,test_ext)



year <- data.frame(water_year=factor("2015",levels=levels(snow_output_dat$water_year)))
test_pred2015 <- raster::predict(test_env,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                     filename="~/GIS/snow_test_lme_2015.tif",overwrite=TRUE)
year <- data.frame(water_year=factor("2011",levels=levels(snow_output_dat$water_year)))
test_pred2011 <- raster::predict(test_env,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                             filename="~/GIS/snow_test_lme_2011.tif",overwrite=TRUE)

brks <- seq(0,250,by=50)
nb <- length(brks)-1
cols <- rev(terrain.colors(nb))

par(mfrow=c(1,2))
plot(test_pred2015,breaks=brks,col=cols)
plot(test_pred2011,breaks=brks,col=cols)

##Raster predictions at 3m resolution for 2009 - 2015.
year <- data.frame(water_year=factor("2009",levels=levels(snow_output_dat$water_year)))
lme_pred_2009 <- raster::predict(env_pred,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                                 filename="~/GIS/snow_lme_pred_2009.tif",overwrite=TRUE)
year <- data.frame(water_year=factor("2010",levels=levels(snow_output_dat$water_year)))
lme_pred_2010 <- raster::predict(env_pred,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                                 filename="~/GIS/snow_lme_pred_2010.tif",overwrite=TRUE)
year <- data.frame(water_year=factor("2011",levels=levels(snow_output_dat$water_year)))
lme_pred_2011 <- raster::predict(env_pred,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                                 filename="~/GIS/snow_lme_pred_2011.tif",overwrite=TRUE)
year <- data.frame(water_year=factor("2012",levels=levels(snow_output_dat$water_year)))
lme_pred_2012 <- raster::predict(env_pred,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                                 filename="~/GIS/snow_lme_pred_2012.tif",overwrite=TRUE)
year <- data.frame(water_year=factor("2013",levels=levels(snow_output_dat$water_year)))
lme_pred_2013 <- raster::predict(env_pred,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                                 filename="~/GIS/snow_lme_pred_2013.tif",overwrite=TRUE)
year <- data.frame(water_year=factor("2014",levels=levels(snow_output_dat$water_year)))
lme_pred_2014 <- raster::predict(env_pred,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                                 filename="~/GIS/snow_lme_pred_2014.tif",overwrite=TRUE)
year <- data.frame(water_year=factor("2015",levels=levels(snow_output_dat$water_year)))
lme_pred_2015 <- raster::predict(env_pred,model=snowmod_lme8,const=year,progress="text",re.form=~0,
                                  filename="~/GIS/snow_lme_pred_2015.tif",overwrite=TRUE)


####Constructs raster predictions for each year based on ordinary least-squares####

##Best linear model
snowmod2013 <- lm(snow_diss_DOY~elevation*water_year + I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m,
                  data=snow_output_dat)
snowmod2012 <- lm(snow_diss_DOY~elevation*water_year + I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m,
                  data=snow_output_dat)
snowmod2011 <- lm(snow_diss_DOY~elevation*water_year + I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m,
                  data=snow_output_dat)
snowmod2010 <- lm(snow_diss_DOY~elevation*water_year + I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m,
                  data=snow_output_dat)
snowmod2009 <- lm(snow_diss_DOY~elevation*water_year + I(elevation^2)+X_UTM*Y_UTM+canopy_pct+slope+relev_30m,
                  data=snow_output_dat)

##Tests the model on a small extent.
test_ext <- extent(c(env_pred@extent@xmin+11000,
                     env_pred@extent@xmin+12000,
                     env_pred@extent@ymin+6000,
                     env_pred@extent@ymin+12000))

##Function to extract estimates and prediction standard errors
predfun <- function(model, data) {
  v <- predict(model, data,interval="prediction",level=0.95)
  se <- v[,3] - v[,2]
  cbind(p=v[,1],lwr=v[,2],upr=v[,3])
}

##Function to extract estimates and prediction standard errors
predfun2 <- function(model, data) {
  v <- predict(model, data,interval="prediction",level=0.50)
  se <- v[,3] - v[,2]
  cbind(p=v[,1],lwr=v[,2],upr=v[,3])
}

##Function to extract estimates and ordinary standard errors
predfun3 <- function(model, data) {
  v <- predict(model, data,se.fit=TRUE)
}

test_pred <- predict(env_pred,snowmod2013,
                     progress='text', 
                     ext=test_ext,
                     fun=predfun3,index=1:2,
                     filename="~/GIS/MORA_snow_pred_test.tif",
                     overwrite=T)

## Predicts grids based on the models.
rasterOptions(maxmemory=5e+08,chunksize=1e+07,timer=TRUE)
env_pred$water_year <- 2011
snow_pred_2013 <- predict(env_pred, model=snowmod2013, 
                          fun=predfun, 
                          index=1:3,
                          progress='text',
                          filename="~/GIS/MORA_snow_pred_2013_3m.tif",
                          overwrite=T,
                          datatype="INT2S")
snow_pred_2013_se <- predict(env_pred, model=snowmod2013, 
                          fun=predfun3, 
                          index=1:2,
                          progress='text',
                          filename="~/GIS/MORA_snow_pred_2013_se_3m.tif",
                          overwrite=T,
                          datatype="INT2S")
env_pred$water_year <- 2012                                           
snow_pred_2012 <- predict(env_pred, model=snowmod2012, 
                          fun=predfun, 
                          index=1:3,
                          progress='text',
                          filename="~/GIS/MORA_snow_pred_2012_3m.tif",
                          overwrite=T,
                          datatype="INT2S")
env_pred$water_year <- 2011
snow_pred_2011 <- predict(env_pred, model=snowmod2011, 
                          fun=predfun, 
                          index=1:3,
                          progress='text',
                          filename="~/GIS/MORA_snow_pred_2011_3m.tif",
                          overwrite=T,
                          datatype="INT2S")
env_pred$water_year <- 2010
snow_pred_2010 <- predict(env_pred, model=snowmod2010, 
                          fun=predfun, 
                          index=1:3,
                          progress='text',
                          filename="~/GIS/MORA_snow_pred_2010_3m.tif",
                          overwrite=T,
                          datatype="INT2S")
env_pred$water_year <- 2009
snow_pred_2009 <- predict(env_pred, model=snowmod2009, 
                          fun=predfun, 
                          index=1:3,
                          progress='text',
                          filename="~/GIS/MORA_snow_pred_2009_3m.tif",
                          overwrite=T,
                          datatype="INT2S")

####Builds models that incorporate the Snow-17 forecast
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Snow/statmodel/")
library(R.matlab)

##brings in the forecast
lowres_dem <- raster("./data/usgs_500m_dem.tif")
temp <- raster("./data/forecast1.tif")
s17 <- readMat("./data/S17_SDD_WY1990_WY2013_CALopt1.mat")
s17.doy <- brick(s17$SDD.doy,crs=temp@crs)
s17.doy@extent <- temp@extent
names(s17.doy) <- paste("snow17_sdd_",1990:2013,sep="")

##subsets layers for 2009 - 2013.
s17.yr <- s17.doy[[c(20:24)]] 
plot(s17.yr[["snow17_sdd_2013"]])
s17.yr[s17.yr > 360] <- NA

##Puts the snow data in a spatial points data frame
snow_sp <- snow_output_dat
coordinates(snow_sp) <- ~X_UTM+Y_UTM
proj4string(snow_sp) <- dem@crs
points(snow_sp,add=TRUE)
snow_sp@data$X_UTM <- coordinates(snow_sp)[,1] 
snow_sp@data$Y_UTM <- coordinates(snow_sp)[,2]


##Samples the snow17 forecast at the snow sensor locations.
snow_samp <- extract(s17.yr,snow_sp,sp=TRUE,method="bilinear")
snow_samp2 <- extract(lowres_dem,snow_samp,method="bilinear",sp=TRUE)
snow_samp2@data$dem_diff <- snow_samp2@data$usgs_500m_dem - snow_samp@data$elevation

##Gets the forecasts into the appropriate columns for each year
snow_samp2@data$s17_sdd <- NA
snow_samp2@data$s17_sdd[snow_samp2@data$water_year==2013] <- snow_samp2@data$snow17_sdd_2013[snow_samp2@data$water_year==2013]
snow_samp2@data$s17_sdd[snow_samp2@data$water_year==2012] <- snow_samp2@data$snow17_sdd_2012[snow_samp2@data$water_year==2012]
snow_samp2@data$s17_sdd[snow_samp2@data$water_year==2011] <- snow_samp2@data$snow17_sdd_2011[snow_samp2@data$water_year==2011]
snow_samp2@data$s17_sdd[snow_samp2@data$water_year==2010] <- snow_samp2@data$snow17_sdd_2010[snow_samp2@data$water_year==2010]
snow_samp2@data$s17_sdd[snow_samp2@data$water_year==2009] <- snow_samp2@data$snow17_sdd_2009[snow_samp2@data$water_year==2009]

##Quick plot to make sure we are on the right track.
##Quick plots to check data.
qplot(snow_diss_DOY,s17_sdd,data=snow_samp2@data,color=Y_UTM)+
  facet_wrap(facets=~water_year)+
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y~poly(x,2))+
  geom_abline(slope=1,intercept=0,linetype=2)+
  labs(x="Snow Dissapearance Date (Julian Day)",
       y="Snow 17 Prediction",
       color="UTM \n North (m)",
       title="Water Year")+
  theme_bw()

snowmod_s17 <- lm(snow_diss_DOY~s17_sdd,data=snow_samp2@data)
summary(snowmod_s17)

summary(snow17_stat1)
