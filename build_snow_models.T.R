##Script to build statistical models of snow duration.
##Author:Ian Breckheimer

####Sets up the workspace####

## Loads required packages
library(lme4)
library(ggmap)
library(raster)
library(ggplot2)
library(geoR)

## Sets working directory
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/cleaned")

##Loads snow duration data.
snow_output <- read.csv("microclimate_snow_coords_final_2_12_14.csv")
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
qplot(snow_diss_DOY,elevation,data=snow_output,facets=water_year~.,color=canopy_pct*100)+
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y~poly(x,2))+
  labs(x="Snow Dissapearance Date (Julian Day)",
       y="Elevation (m)",
       color="Canopy\nCover (%)",
       title="Water Year")+
  theme_bw()

####Preliminary statistical models of SDD####

##Only 2 obs for 2007 so remove.
snow_output_dat <- subset(snow_output,water_year %in% c("2009","2010","2011","2012","2013"))

####Replot with 2007 removed.
##Quick plots to check data.

##Color scale consistent with map.
snow_cols <- c(rgb(250/255,8/255,8/255),
               rgb(241/255,165/255,191/255),
               rgb(123/255,172/255,30/255),
               rgb(50/255,98/255,132/255),
               rgb(136/255,7/255,129/255),
               rgb(0,0,0))

# ##Melt dates for paradise SNOTEL 2009 - 2013
# para_dates <- c(202,205,240,209,199)
# para_elevs <- rep(1563,5)
# 
# ##Appends those values to data frame.
# snow_output_para <- snow_output_dat[1:5,]
# snow_output_para[,] <- rep(NA,25)
# snow_output_para$snow_diss_DOY <- para_dates
# snow_output_para$elevation <- para_elevs
# snow_output_para$study <- rep("Para_SNOTEL",5)
# snow_output_para$site_name <- rep("Para_SNOTEL",5)
# snow_output_para$water_year <- 2009:2013
# snow_output_para$PARA <- TRUE

# snow_output_dat$PARA <- FALSE
# snow_output_dat <- rbind(snow_output_dat,snow_output_para)


cairo_pdf("sdd_plot_allyears.pdf",width=2.5,height=5,family="CMU Serif")
ggplot(data=snow_output_dat)+
  facet_grid(water_year~.)+
  geom_point(aes(x=snow_diss_DOY,y=elevation,color=snow_diss_DOY,size=PARA,shape=PARA))+
  geom_smooth(aes(x=snow_diss_DOY,y=elevation),method = "lm", se=TRUE, color="black", formula = y~poly(x,2))+
  scale_colour_gradientn(colours=snow_cols, values = NULL,
                         space = "Lab", na.value = "grey50",
                         guide = "colourbar")+
  scale_size_manual(values=c(1,4))+
  scale_x_continuous(breaks=c(91,152,213), 
                     labels=c("April", "June","August"))+
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
snowmod9 <- lm(snow_diss_DOY~elevation*water_year+I(elevation^2)+X_UTM*Y_UTM+elevation*relev_30m,data=snow_output_dat)
summary(snowmod9)
AIC(snowmod,snowmod2,snowmod3,snowmod4,snowmod5,snowmod6,snowmod7,snowmod8,snowmod9)

####Examines spatial autocorrelation in the residuals####

## Subsets the data
snow_spat <- subset(snow_output,water_year %in% c("2009","2010","2011","2012","2013"))

## Residuals from best model
snow_spat$resid <- snowmod8$residuals

## Converts coordinates to km to avoid numerical problems with fitting variograms
obs_coords <- as.matrix(snow_spat[,c("X_UTM","Y_UTM")])/1000

## Converts data to a spatial points data frame.
coordinates(snow_spat) <- obs_coords

## Computes the empirical variogram on the raw data and residuals, 
## then fits a variogram model to the data.
snow_vg <- variog(coords=snow_spat@coords,data=snow_spat@data$snow_diss_DOY,option="bin",
                  breaks=seq(0,8,by=0.01),max.dist=2)
snow_fit <- variofit(snow_vg,ini.cov.pars=c(400,0.2),nugget=0,cov.model="exp",fix.nugget=FALSE)

resid_vg <- variog(coords=snow_spat@coords,data=snow_spat@data$resid,option="bin",
                   breaks=seq(0,8,by=0.01),max.dist=2)
resid_fit <- variofit(resid_vg,ini.cov.pars=c(200,0.6),nugget=0,cov.model="exp",fix.nugget=FALSE)

## Plots the variogram models.
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/maps")
svg("snow_variogram.svg",width=4,height=4,family='sans')
par(mfrow=c(2,1),mar=c(0,3,0,0),oma=c(4,1,1,1))
plot(snow_vg,ylim=c(0,480),pch=19,col=rgb(0,0.7,0.9,0.5),cex=0.5,xaxt='n')
lines.variomodel(snow_fit,lwd=1.5)
text(2,50,labels="Raw Data",pos=2)
plot(resid_vg,pch=19,col=rgb(0,0.7,0.2,0.5),cex=0.5,ylim=c(0,480))
lines.variomodel(resid_fit,lwd=1.5)
text(2,400,labels="Model Residuals",pos=2)
mtext("Distance (km)",side=1,outer=T,padj = 4)
mtext("Semivariance",side=2,outer=T,padj = 1)
dev.off()

####Mixed-effects models to deal with non-independence by location####
snowmod_lme <- lmer(snow_diss_DOY~elevation + water_year + (1 |site_name),
                    data=snow_output_dat,REML=FALSE)
summary(snowmod_lme)
snowmod_lme2 <- lmer(snow_diss_DOY~elevation + (1 + water_year|site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme2)
snowmod_lme3 <- lmer(snow_diss_DOY~elevation + canopy_pct + (1 + water_year|site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme3)
snowmod_lme4 <- lmer(snow_diss_DOY~elevation + canopy_pct + Y_UTM + (1 + water_year|site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme4)
snowmod_lme5 <- lmer(snow_diss_DOY~elevation + canopy_pct + relev_30m + Y_UTM + (1 + water_year|site_name),
                     data=snow_output_dat,REML=FALSE)
summary(snowmod_lme5)
snowmod_lme6 <- lmer(snow_diss_DOY~water_year+elevation + I(elevation^2) + canopy_pct + slope + relev_30m + Y_UTM + X_UTM + (1 + water_year|site_name),
                     data=snow_output_dat,REML=FALSE,control=lmerControl(optCtrl=list(maxfun=20000)))
summary(snowmod_lme6)
snowmod_lme7 <- lmer(snow_diss_DOY~elevation*water_year + I(elevation^2) + slope + canopy_pct + relev_30m*elevation + Y_UTM + X_UTM + (1 + water_year|site_name),
                     data=snow_output_dat,REML=FALSE,control=lmerControl(optCtrl=list(maxfun=20000)))
summary(snowmod_lme7)
resid <- residuals(snowmod_lme7)
plot(density(resid))
AIC(snowmod_lme,snowmod_lme2,snowmod_lme3,snowmod_lme4,snowmod_lme5,snowmod_lme6,snowmod_lme7)

##Best lme model with REML for estimating parameters
snowmod_lme8 <- lmer(snow_diss_DOY~water_year*elevation + I(elevation^2) + slope + canopy_pct + relev_30m*elevation + X_UTM + Y_UTM + (1 + water_year|site_name),
                     data=snow_output_dat,REML=TRUE,control=lmerControl(optCtrl=list(maxfun=20000)))
summary(snowmod_lme8)

##Puts the lme residuals on a map
snow_output_dat$lmeresid <- residuals(snowmod_lme8)
map <- get_map(location=c(-121.75,46.85),zoom=11,maptype="satellite",source="google")
ggmap(map)+
  geom_point(aes(x = jitter(longitude,amount=0.005), 
                 y = jitter(latitude,amount=0.005), 
                 color = lmeresid), 
             data = subset(snow_output_dat,lmeresid > 10 | lmeresid < -10) , alpha = .8)+
  scale_colour_gradientn(colours = c("#EDAAAA","#FFFFFF","#0800FF"),values=c(0,0.5,1))

# env_pred <- brick(list(elevation=env3$elev,canopy_pct=env3$can_pct,slope=env3$slope))
# gf <- focalWeight(env_pred$elevation, 30, "circle")
# env_pred$relev_30m <- env_pred$elevation - focal(env_pred$elevation,gf)
# env_pred$water_year <- 2013
# env_pred$X_UTM <- raster('~/GIS/MORA_UTMX_3m.tif')
# env_pred$Y_UTM <- raster('~/GIS/MORA_UTMY_3m.tif')
# writeRaster(env_pred,"~/GIS/MORA_snow_env_pred2.grd",format="raster",overwrite=T)
env_pred <- brick("~/GIS/MORA_snow_env_pred2.grd")
env_pred <- env_pred[[c(1:4,6:7)]]

###Constructs raster predictions for each year based on the best mixed-effects model
## Subsets the data
snow_output_dat <- subset(snow_output,water_year %in% c("2009","2010","2011","2012","2013"))
pred <- predict(snowmod_lme8,REform = NA)
plot(pred~snow_output_dat$snow_diss_DOY)

pred2013 <- function(terms){
  coefs <- fixef(snowmod_lme8)
  pred <- coefs[1]+coefs[5]+coefs[6]*terms[[1]]+coefs[7]*terms[[1]]^2+
    coefs[8]*terms[[3]]+coefs[9]*terms[[2]]+coefs[10]*terms[[4]]+coefs[11]*terms[[5]]+
    coefs[12]*terms[[6]]+coefs[16]*terms[[1]]+coefs[17]*terms[[1]]*terms[[4]]
  return(pred)
}

pred2012 <- function(terms){
  coefs <- fixef(snowmod_lme8)
  pred <- coefs[1]+coefs[4]+coefs[6]*terms[[1]]+coefs[7]*terms[[1]]^2+
    coefs[8]*terms[[3]]+coefs[9]*terms[[2]]+coefs[10]*terms[[4]]+coefs[11]*terms[[5]]+
    coefs[12]*terms[[6]]+coefs[15]*terms[[1]]+coefs[17]*terms[[1]]*terms[[4]]
  return(pred)
}

pred2011 <- function(terms){
  coefs <- fixef(snowmod_lme8)
  pred <- coefs[1]+coefs[3]+coefs[6]*terms[[1]]+coefs[7]*terms[[1]]^2+
    coefs[8]*terms[[3]]+coefs[9]*terms[[2]]+coefs[10]*terms[[4]]+coefs[11]*terms[[5]]+
    coefs[12]*terms[[6]]+coefs[14]*terms[[1]]+coefs[17]*terms[[1]]*terms[[4]]
  return(pred)
}

pred2010 <- function(terms){
  coefs <- fixef(snowmod_lme8)
  pred <- coefs[1]+coefs[2]+coefs[6]*terms[[1]]+coefs[7]*terms[[1]]^2+
    coefs[8]*terms[[3]]+coefs[9]*terms[[2]]+coefs[10]*terms[[4]]+coefs[11]*terms[[5]]+
    coefs[12]*terms[[6]]+coefs[13]*terms[[1]]+coefs[17]*terms[[1]]*terms[[4]]
  return(pred)
}

pred2009 <- function(terms){
  coefs <- fixef(snowmod_lme8)
  pred <- coefs[1]+coefs[6]*terms[[1]]+coefs[7]*terms[[1]]^2+
    coefs[8]*terms[[3]]+coefs[9]*terms[[2]]+coefs[10]*terms[[4]]+coefs[11]*terms[[5]]+
    coefs[12]*terms[[6]]+coefs[17]*terms[[1]]*terms[[4]]
  return(pred)
}

##Tests the model on a small extent.
test_ext <- extent(c(env_pred@extent@xmin+11000,
                     env_pred@extent@xmin+12000,
                     env_pred@extent@ymin+6000,
                     env_pred@extent@ymin+12000))
test_env <- crop(env_pred,test_ext)
test_pred <- overlay(test_env,fun=pred2013,filename="~/GIS/MORA_snow_pred_lme_test.tif",
                     overwrite=TRUE,progress='text')
plot(test_pred)

##Raster predictions at 3m resolution for 2013.
rasterOptions(maxmemory=5e+08,chunksize=1e+07,timer=TRUE)
# lme_pred_2013 <- overlay(x=env_pred,fun=pred2013,
#                          filename="~/GIS/MORA_snow_pred_2013_lme.tif",
#                          overwrite=TRUE,
#                          progress="text")
# 
# ##Raster predictions at 3m resolution for 2012.
# lme_pred_2012 <- overlay(x=env_pred,fun=pred2012,
#                          filename="~/GIS/MORA_snow_pred_2012_lme.tif",
#                          overwrite=TRUE,
#                          progress="text")
# 
# ##Raster predictions at 3m resolution for 2011.
# lme_pred_2011 <- overlay(x=env_pred,fun=pred2011,
#                          filename="~/GIS/MORA_snow_pred_2011_lme.tif",
#                          overwrite=TRUE,
#                          progress="text")

##Raster predictions at 3m resolution for 2010.
lme_pred_2010 <- overlay(x=env_pred,fun=pred2010,
                         filename="~/GIS/MORA_snow_pred_2010_lme.tif",
                         overwrite=TRUE,
                         progress="text")

##Raster predictions at 3m resolution for 2009.
lme_pred_2009 <- overlay(x=env_pred,fun=pred2009,
                         filename="~/GIS/MORA_snow_pred_2009_lme.tif",
                         overwrite=TRUE,
                         progress="text")

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
