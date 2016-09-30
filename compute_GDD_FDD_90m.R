##Script to calculate GDD, FDD, and snow-free GDD from daily data.
library(raster)
library(doParallel)
library(foreach)
rasterOptions(maxmemory=4e+08 )

airtemp_path <- "/Volumes/ib_working/MORAclim_bayes/"
snow_path <- "/Volumes/ib_working/MORAclim_snow/"

snowdayseq <- seq(as.Date("2009-09-01"),as.Date("2015-09-30"),by="day")
airdayseq <- seq(as.Date("2009-01-01"),as.Date("2015-09-30"),by="day")
measdayseq <- seq(as.Date("2009-09-01"),as.Date("2015-08-30"),by="day")
gdd_doy <- as.numeric(strftime(airdayseq,format="%j"))
gdd_yn <- gdd_doy >= 91 & gdd_doy <= 273

##Determines which files to include in which analyses.

gdd_days <- paste("day_",91:273,sep="")
files <- list.files(airtemp_path,pattern="_est.tif$")
tavg_files <- grep("tavg",files,value=TRUE)
tavg_gdd <- tavg_files[gdd_yn]

tmin_files <- grep("tmin",files,value=TRUE)
tmin_gdd <- tmin_files[gdd_yn]

tmax_files <- grep("tmax",files,value=TRUE)
tmax_gdd <- tmax_files[gdd_yn]

snow_files <- list.files(snow_path,pattern=".tif$")
snow_dates <- sub(".tif","",sub("snowprob_","",snow_files))
snow_doys <- as.numeric(strftime(as.Date(snow_dates),format="%j"))
snow_yn <- snow_doys >= 91 & snow_doys <= 273
snow_gdd <- snow_files[snow_yn]

####Calculates accumulated GDD,FDD,and snow-free GDD for each year, and averages across years.####
gdd_fun <- function(x){
  x[x < 50000] <- 0
  sum(x/10000)
}
sfgdd_fun <- function(x){
  x[x < 5] <- 0
  sum(x)
}
fdd_fun <- function(x){
  x[x > 0] <- 0
  sum(x/10000)*-1
}

sf_fun <- function(x,y){
  return((x/10000)*((y/10000)*-1+1))
}

gmin_fun <- function(x){
  return(min(x/10000))
}

gmax_fun <- function(x){
  return(max(x/10000))
}

ffd_fun <- function(x){
  ifelse(any(x<0),min(which(x<0)),NA) + 212
}

lfd_fun <- function(x){
  ifelse(any(x<0),max(which(x<0)),NA)
}


setwd(airtemp_path)
gdd_temp <- raster(tavg_gdd[1])
gdd_temp[] <- NA
gdd_sum <- stack(gdd_temp)
fdd_sum <- stack(gdd_temp)
sf_gdd_sum <- stack(gdd_temp)
sf_fdd_sum <- stack(gdd_temp)
gmin <- stack(gdd_temp)
gmax <- stack(gdd_temp)
ffd <- stack(gdd_temp)
lfd <- stack(gdd_temp)

years <- 2009:2015

for(i in 1:length(years)){
  print(paste("Processing data from year", years[i]))
  year_files <- grep(as.character(years[i]),tavg_gdd,fixed=TRUE,value=TRUE)
  year_snow_files <- grep(as.character(years[i]),snow_gdd,fixed=TRUE,value=TRUE)
  year_tmin_files <-  grep(as.character(years[i]),tmin_gdd,fixed=TRUE,value=TRUE)
  year_tmin_files_all <-  grep(as.character(years[i]),tmin_files,fixed=TRUE,value=TRUE)
  year_tmin_files_spr <- year_tmin_files_all[1:212]
  year_tmin_files_fall <- year_tmin_files_all[213:length(year_tmin_files_all)]
  year_tmax_files <-  grep(as.character(years[i]),tmax_gdd,fixed=TRUE,value=TRUE)

  # year_stack <- stack(year_files,quick=TRUE)
  # year_brick <- brick(year_stack)
  # 
  # year_tmin_stack <- stack(year_tmin_files,quick=TRUE)
  # year_tmin_brick <- brick(year_tmin_stack)
  
  spring_tmin_stack <- stack(year_tmin_files_spr,quick=TRUE)
  spring_tmin_brick <- brick(spring_tmin_stack)
  fall_tmin_stack <- stack(year_tmin_files_fall,quick=TRUE)
  fall_tmin_brick <- brick(fall_tmin_stack)
  
  # year_tmax_stack <- stack(year_tmax_files,quick=TRUE)
  # year_tmax_brick <- brick(year_tmax_stack)
  # 
  # snow_stack <- stack(paste("../MORAclim_snow/",year_snow_files,sep=""),quick=TRUE)
  # snow_brick <- brick(snow_stack)
  # 
  # gdd_sum[[i]] <- calc(year_brick,fun=gdd_fun,progress="text")
  # fdd_sum[[i]] <- calc(year_tmin_brick,fun=fdd_fun,progress="text")
  # sf_temp <- overlay(year_brick,snow_brick,fun=sf_fun,progress="text")
  # sf_tmin <- overlay(year_tmin_brick,snow_brick,fun=sf_fun,progress="text")
  # sf_gdd_sum[[i]] <- calc(sf_temp,fun=sfgdd_fun,progress="text")
  # sf_fdd_sum [[i]] <- calc(sf_tmin,fun=fdd_fun,progress="text")
  # gmin[[i]] <- calc(year_tmin_brick,fun=gmin_fun)
  # gmax[[i]] <- calc(year_tmax_brick,fun=gmax_fun)
  ffd[[i]] <- calc(fall_tmin_brick,fun=ffd_fun)
  lfd[[i]] <- calc(spring_tmin_brick,fun=lfd_fun)
}
# names(gdd_sum) <- paste("Year",years)
# names(fdd_sum) <- paste("Year",years)
# names(sf_gdd_sum) <- paste("Year",years)
# names(sf_fdd_sum) <- paste("Year",years)
# names(gmin) <- paste("Year",years)
# names(gmax) <- paste("Year",years)
names(ffd) <- paste("Year",years)
names(lfd) <- paste("Year",years)

##Computes yearly growing-season averages.
# gdd_avg <- calc(gdd_sum,fun=function(x){mean(x)},overwrite=TRUE,
#                 filename="../MORAclim_summaries/gdd_avg_2010_2015.tif")
# fdd_avg <- calc(fdd_sum,fun=function(x){mean(x)},overwrite=TRUE,
#                 filename="../MORAclim_summaries/fdd_avg_2010_2015.tif")
# sf_gdd_avg <- calc(sf_gdd_sum,fun=function(x){mean(x)},overwrite=TRUE,
#                    filename="../MORAclim_summaries/sfgdd_avg_2010_2015.tif")
# sf_fdd_avg <- calc(sf_fdd_sum,fun=function(x){mean(x)},overwrite=TRUE,
#                    filename="../MORAclim_summaries/sffdd_avg_2010_2015.tif")
# gmin_avg <- calc(gmin,fun=function(x){mean(x)},overwrite=TRUE,
#                  filename="../MORAclim_summaries/gsmin_avg_2010_2015.tif")
# gmax_avg <- calc(gmax,fun=function(x){mean(x)},overwrite=TRUE,
#                  filename="../MORAclim_summaries/gsmax_avg_2010_2015.tif")
lfd_avg <- calc(lfd,fun=function(x){mean(x)},overwrite=TRUE,
                 filename="../MORAclim_summaries/lfd_avg_2009_2015.tif")
ffd_avg <- calc(ffd[[1:6]],fun=function(x){mean(x)},overwrite=TRUE,
                 filename="../MORAclim_summaries/ffd_avg_2009_2014.tif")

plot(gdd_sum,zlim=c(0,3000))
plot(fdd_sum,zlim=c(0,2000))
plot(sf_gdd_sum,zlim=c(0,3000))
plot(sf_fdd_sum,zlim=c(0,50))

##Calculates snow cover stats for each water year.
wy_days <- seq(as.Date("2009-09-01"),as.Date("2015-09-22"),by="day")
wyears <- rep(NA,length(wy_days))
wyears[wy_days %in% seq(as.Date("2009-09-15"),as.Date("2010-09-14"),by='day')] <- 2010
wyears[wy_days %in% seq(as.Date("2010-09-15"),as.Date("2011-09-14"),by='day')] <- 2011
wyears[wy_days %in% seq(as.Date("2011-09-15"),as.Date("2012-09-14"),by='day')] <- 2012
wyears[wy_days %in% seq(as.Date("2012-09-15"),as.Date("2013-09-14"),by='day')] <- 2013
wyears[wy_days %in% seq(as.Date("2013-09-15"),as.Date("2014-09-14"),by='day')] <- 2014
wyears[wy_days %in% seq(as.Date("2014-09-15"),as.Date("2015-09-14"),by='day')] <- 2015

setwd(snow_path)
snow_temp <- raster(snow_files[1])
snow_temp[] <- NA
snow_days <- stack(gdd_temp)
first_snow <- stack(gdd_temp)
last_snow <- stack(gdd_temp)

last_snow_fun <- function(x){
  max(which(x>5000))
}

first_snow_fun <- function(x){
  min(which(x>5000))
}

wy <- 2010:2015
for(i in 1:length(wy)){
  print(paste("Processing Water Year ",wy[i],sep=""))
  snow_files_wy <- snow_files[which(wyears == wy[i])]
  snow_stack_wy <- stack(paste(snow_path,snow_files_wy,sep=""),quick=TRUE)
  snow_brick_wy <- brick(snow_stack_wy)
  
  ##Removes problem day in 2015
  if(wy[i]==2015){
    snow_brick_wy[[356]] <- 0 
  }
  
  snow_days[[i]] <- calc(snow_brick_wy,fun=sum) / 10000
  first_snow[[i]] <- calc(snow_brick_wy,fun=first_snow_fun)
  last_snow[[i]] <- calc(snow_brick_wy,fun=last_snow_fun)
}

##Calculates 2010 - 2015 averages.
scd_mean <- calc(snow_days,fun=mean,
                 filename="../MORAclim_summaries/scd_avg.tif",
                 overwrite=TRUE)
lsd_mean <- calc(last_snow,fun=mean)
fsd_mean <- calc(first_snow,fun=mean)

##Puts lsd and fsd back in terms of day of year.
fsd_mean_doy <- fsd_mean + 258
writeRaster(fsd_mean_doy,filename="../MORAclim_summaries/fsd_avg_2009_2014.tif",
            overwrite=TRUE)
lsd_mean_doy <- lsd_mean - 107
writeRaster(lsd_mean_doy,filename="../MORAclim_summaries/lsd_avg_2010_2015.tif")

####Calculates seasonal and annual temperature averages.
files <- list.files(airtemp_path,pattern="_est.tif$")
tmax_files <- grep("tmax",files,fixed=TRUE,value=TRUE)
tmin_files <- grep("tmin",files,fixed=TRUE,value=TRUE)
tavg_files <- grep("tavg",files,fixed=TRUE,value=TRUE)
process_yearly_yn <- airdayseq %in% measdayseq

seas_tmax_files <- tmax_files[process_yearly_yn]
seas_tmin_files <- tmin_files[process_yearly_yn]
seas_tavg_files <- tavg_files[process_yearly_yn]

measmonthseq <- as.numeric(strftime(measdayseq,format="%m"))
djf_days <- measmonthseq %in% c(12,1,2)
mam_days <- measmonthseq %in% c(3,4,5)
jja_days <- measmonthseq %in% c(6,7,8)
son_days <- measmonthseq %in% c(9,10,11)

djf_tmax_files <- seas_tmax_files[djf_days]
mam_tmax_files <- seas_tmax_files[mam_days]
jja_tmax_files <- seas_tmax_files[jja_days]
son_tmax_files <- seas_tmax_files[son_days]

djf_tmin_files <- seas_tmin_files[djf_days]
mam_tmin_files <- seas_tmin_files[mam_days]
jja_tmin_files <- seas_tmin_files[jja_days]
son_tmin_files <- seas_tmin_files[son_days]

djf_tavg_files <- seas_tavg_files[djf_days]
mam_tavg_files <- seas_tavg_files[mam_days]
jja_tavg_files <- seas_tavg_files[jja_days]
son_tavg_files <- seas_tavg_files[son_days]

avg_fun <- function(x){mean(x)/10000}
sd_fun <- function(x){sd(x)/10000}
quantfun <- function(x,...){quantile(x/10000,probs=c(0.05,0.25,0.75,0.95),na.rm=TRUE)}

compute_seasonal_averages <- function(files,outname,overwrite=FALSE){
  library(raster)
  rasterOptions(maxmemory=16e+08 )
  seas_stack <- stack(files,quick=TRUE)
  seas_brick <- brick(seas_stack)
  print(paste("Calculating averages"))
  avg <- calc(seas_brick,fun=avg_fun)
  sd <- calc(seas_brick,fun=sd_fun)
  print(paste("Calculating quantiles"))
  quants <- calc(seas_brick,fun=quantfun)
  outbrick <- brick(avg,quants,sd,filename=outname,overwrite=overwrite)
  names(outbrick) <- c("avg","qt05","qt25","qt75","qt95","sd")
  return(outbrick)
}

file_list <- list(djf_tmax_files,mam_tmax_files,jja_tmax_files,son_tmax_files,
                  djf_tmin_files,mam_tmin_files,jja_tmin_files,son_tmin_files,
                  djf_tavg_files,mam_tavg_files,jja_tavg_files,son_tavg_files)
outnames <- c("djf_tmax_sum.tif","mam_tmax_sum.tif","jja_tmax_sum.tif","son_tmax_sum.tif",
              "djf_tmin_sum.tif","mam_tmin_sum.tif","jja_tmin_sum.tif","son_tmin_sum.tif",
              "djf_tavg_sum.tif","mam_tavg_sum.tif","jja_tavg_sum.tif","son_tavg_sum.tif")
outnames <- paste("../MORAclim_summaries/",outnames,sep="")

setwd("/Volumes/ib_working/MORAclim_bayes/")
cl <- makeCluster(3)
registerDoParallel(cl)
summaries <- foreach(i = 1:length(outnames)) %dopar% {
  compute_seasonal_averages(file_list[[i]],outname=outnames[i],overwrite=TRUE)
}
stopCluster(cl)

##Computes annual averages and quantiles.
ann_file_list <- list(seas_tmax_files,seas_tmin_files,seas_tavg_files)
ann_outnames <- c("ann_tmax_sum.tif","ann_tmin_sum.tif","ann_tavg_sum.tif")
ann_outnames <- paste("../MORAclim_summaries/",ann_outnames,sep="")

setwd("/Volumes/ib_working/MORAclim_bayes/")
cl <- makeCluster(3)
registerDoParallel(cl)
ann_summaries <- foreach(i = 1:length(ann_outnames)) %do% {
  compute_seasonal_averages(ann_file_list[[i]],outname=ann_outnames[i],overwrite=TRUE)
}
stopCluster(cl)

##Computes average snow duration and last snow date


##Explores geographic relationships to see if everything looks reasonable.

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


plot(preds$elev,summaries[[1]][[1]],xlim=c(-200,4500),xlab="",
     ylab="Temperature (C)",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tmax")
plot(preds$elev,summaries[[1]][[2]],xlim=c(-200,4500),xlab="",
     ylab="Temperature (C)",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tmax")
plot(preds$elev,summaries[[1]][[3]],xlim=c(-200,4500),xlab="",add=TRUE,
     ylab="Temperature (C)",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tmax")
plot(preds$elev,summaries[[1]][[4]],xlim=c(-200,4500),xlab="",
     ylab="Temperature (C)",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tmax")
plot(preds$elev,summaries[[1]][[5]],xlim=c(-200,4500),xlab="",add=TRUE,
     ylab="Temperature (C)",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tmax")


##Plots lapse rates by season.
pdf("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/MORAclim_lapse_seas.pdf",
    width=6,height=4)
par(mfrow=c(1,3),mar=c(2,2,0,0))
##Tmax
plot(preds$elev,summaries[[1]][[1]],xlim=c(-200,4500),xlab="",
     ylab="Temperature (C)",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tmax")
lm_djf_tmax <- lm(summaries[[1]][[1]][]~preds$elev[])
abline(coefficients(lm_djf_tmax)[1],coefficients(lm_djf_tmax)[2],lty=2,col=rgb(0,0,0))
text(400,0,labels=paste(round(coefficients(lm_djf_tmax)[2]*1000,digits=2),"\n C/km"),
     cex=0.9)

plot(preds$elev,summaries[[2]][[1]],xlim=c(0,4500),add=TRUE,
     col=rgb(0,1,0,0.005),ylim=c(-10,25),main="Tmax")
lm_mam_tmax <- lm(summaries[[2]][[1]][]~preds$elev[])
abline(coefficients(lm_mam_tmax)[1],coefficients(lm_mam_tmax)[2],lty=2,col=rgb(0,0.5,0))
text(400,8,labels=paste(round(coefficients(lm_mam_tmax)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0,0.5,0))

plot(preds$elev,summaries[[3]][[1]],xlim=c(-200,4500),add=TRUE,
     col=rgb(1,0,0,0.005),ylim=c(-10,25),main="Tmax")
lm_jja_tmax <- lm(summaries[[3]][[1]][]~preds$elev[])
abline(coefficients(lm_jja_tmax)[1],coefficients(lm_jja_tmax)[2],lty=2,col=rgb(0.5,0,0))
text(400,25,labels=paste(round(coefficients(lm_jja_tmax)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0.5,0,0))

plot(preds$elev,summaries[[4]][[1]],xlim=c(-200,4500),add=TRUE,
     col=rgb(1,0,1,0.005),ylim=c(-10,25),main="Tmax")
lm_son_tmax <- lm(summaries[[4]][[1]][]~preds$elev[])
abline(coefficients(lm_son_tmax)[1],coefficients(lm_son_tmax)[2],lty=2,col=rgb(0.5,0,0.5))
text(400,15,labels=paste(round(coefficients(lm_son_tmax)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0.5,0,0))

##Tavg
plot(preds$elev,summaries[[9]][[1]],xlim=c(-200,4500),xlab="Elevation (m)",
     ylab="",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tavg")
lm_djf_tavg <- lm(summaries[[9]][[1]][]~preds$elev[])
abline(coefficients(lm_djf_tavg)[1],coefficients(lm_djf_tavg)[2],lty=2,col=rgb(0,0,0))
text(400,-2,labels=paste(round(coefficients(lm_djf_tavg)[2]*1000,digits=2),"\n C/km"),
     cex=0.9)

plot(preds$elev,summaries[[10]][[1]],xlim=c(0,4500),add=TRUE,
     col=rgb(0,1,0,0.005),ylim=c(-10,25),main="")
lm_mam_tavg <- lm(summaries[[10]][[1]][]~preds$elev[])
abline(coefficients(lm_mam_tavg)[1],coefficients(lm_mam_tavg)[2],lty=2,col=rgb(0,0.5,0))
text(400,4,labels=paste(round(coefficients(lm_mam_tavg)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0,0.5,0))

plot(preds$elev,summaries[[11]][[1]],xlim=c(-200,4500),add=TRUE,
     col=rgb(1,0,0,0.005),ylim=c(-10,25),main="")
lm_jja_tavg <- lm(summaries[[11]][[1]][]~preds$elev[])
abline(coefficients(lm_jja_tavg)[1],coefficients(lm_jja_tavg)[2],lty=2,col=rgb(0.5,0,0))
text(400,20,labels=paste(round(coefficients(lm_jja_tavg)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0.5,0,0))

plot(preds$elev,summaries[[12]][[1]],xlim=c(-200,4500),add=TRUE,
     col=rgb(1,0,1,0.005),ylim=c(-10,25),main="Tmax")
lm_son_tavg <- lm(summaries[[12]][[1]][]~preds$elev[])
abline(coefficients(lm_son_tavg)[1],coefficients(lm_son_tavg)[2],lty=2,col=rgb(0.5,0,0.5))
text(400,12,labels=paste(round(coefficients(lm_son_tavg)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0.5,0,0))

##Tmin
plot(preds$elev,summaries[[5]][[1]],xlim=c(-200,4500),xlab="",
     ylab="",col=rgb(0,0,0,0.005),ylim=c(-15,25),main="Tmin")
lm_djf_tmin <- lm(summaries[[5]][[1]][]~preds$elev[])
abline(coefficients(lm_djf_tmin)[1],coefficients(lm_djf_tmin)[2],lty=2,col=rgb(0,0,0))
text(400,-4,labels=paste(round(coefficients(lm_djf_tmin)[2]*1000,digits=2),"\n C/km"),
     cex=0.9)

plot(preds$elev,summaries[[6]][[1]],xlim=c(0,4500),add=TRUE,
     col=rgb(0,1,0,0.005),ylim=c(-10,25),main="")
lm_mam_tmin <- lm(summaries[[6]][[1]][]~preds$elev[])
abline(coefficients(lm_mam_tmin)[1],coefficients(lm_mam_tmin)[2],lty=2,col=rgb(0,0.5,0))
text(400,4,labels=paste(round(coefficients(lm_mam_tmin)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0,0.5,0))

plot(preds$elev,summaries[[7]][[1]],xlim=c(-200,4500),add=TRUE,
     col=rgb(1,0,0,0.005),ylim=c(-10,25),main="")
lm_jja_tmin <- lm(summaries[[7]][[1]][]~preds$elev[])
abline(coefficients(lm_jja_tmin)[1],coefficients(lm_jja_tmin)[2],lty=2,col=rgb(0.5,0,0))
text(400,16,labels=paste(round(coefficients(lm_jja_tmin)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0.5,0,0))

plot(preds$elev,summaries[[8]][[1]],xlim=c(-200,4500),add=TRUE,
     col=rgb(1,0,1,0.005),ylim=c(-10,25),main="")
lm_son_tmin <- lm(summaries[[8]][[1]][]~preds$elev[])
abline(coefficients(lm_son_tmin)[1],coefficients(lm_son_tmin)[2],lty=2,col=rgb(0.5,0,0.5))
text(400,10,labels=paste(round(coefficients(lm_son_tmin)[2]*1000,digits=2),"\n C/km"),
     cex=0.9,col=rgb(0.5,0,0))

##Legend
legend(2000,24,legend=c("DJF","MAM","JJA","SON"),
       lty=c(2,2,2,2),bty="n",text.col=c(rgb(0,0,0),
                                         rgb(0,1,0),
                                         rgb(1,0,0),
                                         rgb(1,0,1)),
        col=c(rgb(0,0,0),
              rgb(0,0.5,0),
              rgb(0.5,0,0),
              rgb(0.5,0,0.5)))
dev.off()

##Interpolates to 10m resolution for display.
mat <- ann_summaries[[3]][[1]]
mat_9m <- disaggregate(mat,fact=10,method="bilinear",
                       filename="/Volumes/ib_working/MORAclim_summaries/mat_9m.tif")
sum_tmax <- summaries[[3]][[1]]
sum_tmax_9m <- disaggregate(sum_tmax,fact=10,method="bilinear",
                            filename="/Volumes/ib_working/MORAclim_summaries/sum_tmax_9m.tif")
win_tmin <- summaries[[5]][[1]]
win_tmin_9m <- disaggregate(win_tmin,fact=10,method="bilinear",
                            filename="/Volumes/ib_working/MORAclim_summaries/win_tmin_9m.tif")

##Plots residuals after removing elevation trend.
rand <- sampleRandom(stack(mat,win_tmin,sum_tmax,preds),size=60000)
mat_elev_lm <- lm(avg.1~elev,data=data.frame(rand))
mat_elev_pred <- predict(preds,mat_elev_lm)
mat_elev_resid <- mat - mat_elev_pred

sum_tmax_lm <- lm(avg.3~elev,data=data.frame(rand))
sum_tmax_elev_pred <- predict(preds,sum_tmax_lm)
sum_tmax_elev_resid <- sum_tmax - sum_tmax_elev_pred
plot(sum_tmax_elev_resid)

win_tmin_lm <- lm(avg.2~elev,data=data.frame(rand))
win_tmin_elev_pred <- predict(preds,win_tmin_lm)
win_tmin_elev_resid <- win_tmin - win_tmin_elev_pred
plot(win_tmin_elev_resid)
resids <- stack(mat_elev_resid,sum_tmax_elev_resid,win_tmin_elev_resid)
names(resids) <- c("Tavg","Sum. Tmax","Win Tmin")
pal <- colorRampPalette(c("red","white","blue"))

plot(resids,zlim=c(-4,4),
     nc=3)

