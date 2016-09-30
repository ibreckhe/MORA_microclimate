##Script to assemble snow cover data for 2009 to 2015.
library(xts)
library(dplyr)

setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/processed/unflagged/")
snowmeta <- read.csv("../Snow_cover_meta_cleaned_geo_10_16_2015.csv")

raw_files <- list.files(".",pattern=".csv")

snow_fun <- function(ts){
  (min(ts) - max(ts)) < 1 & max(ts) < 2
}

##Loops through each file and estimates whether there was snow on each day.
for(i in 1:length(raw_files)){
  print(paste("Now processing file ",i," of ",length(raw_files)))
  d <- read.csv(raw_files[i])
  d <- d[complete.cases(d),]
  d$DATE <- paste(d$YEAR,d$MONTH,d$DAY,sep="-")
  d$TIME <- paste(d$HOUR,":00",sep="")
  d$DATETIME <- as.POSIXct(paste(d$DATE,d$TIME,sep=" "),tz="Etc/GMT-7")
  ts <- xts(d$TEMP,order.by=d$DATETIME)
  ts_min_daily <- apply.daily(ts,FUN=min)
  ts_max_daily <- apply.daily(ts,FUN=max)
  ts_snow <- apply.daily(ts,FUN=snow_fun)
  df <- data.frame(datetime=as.Date(index(ts_snow)),
                   min_temp_C=ts_min_daily,
                   max_temp_C=ts_max_daily,
                   snow=ts_snow)
  write.csv(df,paste("../daily_soiltemp/",raw_files[i],sep=""),row.names=FALSE)
}

##Constructs a daily dataset with all sensors.
setwd("../daily_soiltemp/")
sites <- unique(snowmeta$site_sensor)
daily_files <- list.files(".",pattern=".csv")
days <- seq(as.Date("2009-09-01"),as.Date("2015-09-30"),by="day")
day_ts <- xts(rep(FALSE,length(days)),order.by=days)

for (i in 1:length(sites)){
  print(paste("Processing site",sites[i]))
  site_files <- unique(snowmeta$out_filename[snowmeta$site_sensor==sites[i]])
  site_daily <- daily_files[daily_files %in% site_files]
  fd <- read.csv(site_daily[1])
  ts_all <- xts(fd$snow,order.by=as.Date(fd$datetime))
  if(length(site_daily)>=2){
    for(j in 2:length(site_daily)){
      fd <- read.csv(site_daily[j])
      ts <- xts(fd$snow,order.by=as.Date(fd$datetime))
      ts_all <- c(ts_all,ts)
    }
  }
  
  names(ts_all) <- sites[i]
  day_ts <- merge(day_ts,ts_all,all=c(TRUE,FALSE))
}

day_df <- data.frame(datetime=index(day_ts),day_ts)
DOY <- as.numeric(format(day_df$datetime,format="%j"))
Year <- as.numeric(format(day_df$datetime,format="%Y"))

##Replaces missing summer measurements with zeros.
for(i in 2009:2015){
  for(j in 2:ncol(day_df)){
    prevtest <- !is.na(day_df[Year==i,j][1])
    nexttest <- any(!is.na(day_df[Year==i & DOY >= 273,j]))
    if(prevtest & nexttest){
      daysnow <- day_df[Year==i,j]
      daysnow[is.na(daysnow)] <- 0
      day_df[Year==i,j] <- daysnow
    }
  }
}


day_df <- day_df[,-2]
write.csv(day_df,"../snow_cover_daily_allsites.csv",row.names=FALSE)

##Loops through each day and fits a logistic regression to predict probability of snow.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/processed/unflagged/")
snowmeta <- read.csv("../microclimate_snow_coords_final_10_19_15.csv")
snowmeta$sensor_name2 <- gsub("-",".",snowmeta$sensor_name,fixed=TRUE)
snowmeta$sensor_name2 <- gsub("1791.STR","X1791.STR",fixed=TRUE,snowmeta$sensor_name2)
snowmeta$sensor_name2 <- gsub("1901.STR","X1901.STR",fixed=TRUE,snowmeta$sensor_name2)

##Reads in rasters of predictors
library(raster)
setwd("/Volumes/ib_working/GIS/")
elev3m <- raster("mcpred_elev_3m.tif")
elev <- aggregate(elev3m,fact=30,fun=mean)
coldair3m <-raster("mcpred_reg_coldair_3m.tif")
coldair <- aggregate(coldair3m,fact=30,fun=mean)
srad3m <- raster("mcpred_srad_MAM_3m.tif")
srad <- aggregate(srad3m,fact=30,fun=mean)
slope3m <- raster("mcpred_slope.tif")
slope <- aggregate(slope3m,fact=30,fun=mean)
cancov3m <- raster("MORA_can_pct_3m.tif")
cancov_raw <- aggregate(cancov3m,fact=30,fun=mean)
utmx3m <- raster("mcpred_utmx_3m.tif")
utmx <- aggregate(utmx3m,fact=30,fun=mean)
utmy3m <- raster("mcpred_utmy_3m.tif")
utmy <- aggregate(utmy3m,fact=30,fun=mean)
preds <- stack(elev,coldair,srad,slope,cancov,utmx,utmy)
names(preds) <- c("elev","coldair","srad","slope","cancov","utmx","utmy")
preds$utmx_s <- preds$utmx / 1000
preds$utmy_s <- preds$utmy / 1000
pred_crs <- crs(elev3m)

snow_site_meta <- data.frame(site=snowmeta$sensor_name2,
                             utmx=snowmeta$X_UTM,
                             utmy=snowmeta$Y_UTM)
snow_site_meta <- unique(snow_site_meta)

##Loops through each day and estimates snow cover
for(i in 1:dim(day_df)[1]){
  day_data <- data.frame(site=names(day_df[i,-1]),
                         doy=DOY[i],
                         snow=as.numeric(day_df[i,-1]))
  
  day_meta <- merge(day_data,snow_site_meta,by="site",all=FALSE)
  day_covars <- data.frame(day_meta,extract(preds,cbind(day_meta$utmx,day_meta$utmy),method="simple"))
  day_covars$utmx_scale <- as.numeric(scale(day_covars$utmx))
  utmx_center <- attributes(scale(day_covars$utmx))$'scaled:center'
  utmx_scale <- attributes(scale(day_covars$utmx))$'scaled:scale'
  day_covars$utmy_scale <- as.numeric(scale(day_covars$utmy))
  utmy_center <- attributes(scale(day_covars$utmy))$'scaled:center'
  utmy_scale <- attributes(scale(day_covars$utmy))$'scaled:scale'  
  
  preds$utmx_scale <- (preds$utmx - utmx_center) / utmx_scale
  preds$utmy_scale <- (preds$utmy - utmy_center) / utmy_scale
  
  ##Checks to see if both snow classes span sufficient range of x and y.
  pres <- day_covars[which(day_covars$snow==1),]
  abs <- day_covars[which(day_covars$snow==0),]
  prev <- nrow(pres)/(nrow(pres)+nrow(abs))
  pres_rangex <- range(pres$utmx)
  pres_rangey <- range(pres$utmy)
  abs_rangex <- range(abs$utmx)
  abs_rangey <- range(abs$utmy)
  
  rangex_test_pres <- (pres_rangex[2] - pres_rangex[1]) > 8000
  rangey_test_pres <- (pres_rangey[2] - pres_rangey[1]) > 8000
  rangex_test_abs <- (abs_rangex[2] - abs_rangex[1]) > 8000
  rangey_test_abs <- (abs_rangey[2] - abs_rangey[1]) > 8000
  
  range_test <- rangex_test_pres & rangey_test_pres & rangex_test_abs & rangey_test_abs
  prev_test <- prev < 0.95 & prev > 0.05
  
  if(range_test & prev_test){
    if(DOY[i]>80 & DOY[i] < 250){
      day_glm <- glm(snow~1,family=binomial(link="logit"),data=day_covars)
      day_step <- step(day_glm,scope=snow~I(elev/1000-3.5)+utmx_scale+utmy_scale+slope+I(srad/1e6),
                       direction="forward",k=4)
      day_pred <- predict(preds,day_step,type="response")
    }else{
      day_glm <- glm(snow~1,family=binomial(link="logit"),data=day_covars)
      day_step <- step(day_glm,scope=snow~I(elev/1000-3.5)+utmx_scale+utmy_scale+slope+cancov,
                       direction="forward",k=4)
      day_pred <- predict(preds,day_step,type="response")
    }
  }else{
    print("Insufficient geographic variability, simplifying model.")
    day_glm <- glm(snow~I(elev/1000-3.5),family=binomial(link="logit"),data=day_covars)
    day_step <- step(day_glm,scope=snow~I(elev/1000-3.5),direction="backward",k=4)
    day_pred <- predict(preds,day_step,type="response")
  }

  plot(day_pred,main=day_df[i,1],zlim=c(0,1))
  filename <- paste("/Volumes/ib_working/MORAclim_snow/snowprob_",as.character(day_df[i,1]),".tif",sep="")
  writeRaster(day_pred*10000,filename=filename,datatype="INT2U",overwrite=TRUE)
}

##Creates graphics and stats for model checking.
setwd("/Volumes/ib_working/MORAclim_snow/")


snowmaps <- list.files(".",pattern=".tif$")
pdf("../MORAclim_plots/snow_maps.pdf",width=5,height=5)
for(i in 1:length(snowmaps)){
  print(paste("plotting map ", i, " of ",length(snowmaps)))
  map <- raster(snowmaps[i]) / 10000
  plot(map,main=snowmaps[i],zlim=c(0,1))
}
dev.off()
