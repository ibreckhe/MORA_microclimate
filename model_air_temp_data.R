## Script to fill gaps in air temperature data.
## Author: Ian Breckheimer
## Date: 16 October 2014

#### Sets up workspace and load data ####

## Loads required packages.
library(raster)
library(spBayes)
library(MBA)
library(geoR)
library(plotrix)
library(fields)
library(boot)

## Sources my functions.
source("~/code/MORA_microclimate/air_temp_functions.R")

## Loads and formats data.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
topowx <- read.csv("topowx_daily_2004_2012.csv")
topowx$DATE <- as.POSIXct(topowx$DATE)
topowx$TAVG <- (topowx$TMAX + topowx$TMIN)/2

fsites <- read.csv("franklindata_utm_covars.csv")
fsites <- fsites[complete.cases(fsites),]
fsites$SITE <- paste(fsites$SITE_NAME,fsites$PLOT_NO,sep="")
fsites <- fsites[order(fsites$SITE),]

prism <- read.csv("PRISM_daily_2004_2015.csv")
prism$DATE <- as.POSIXct(prism$DATE)

topowx_tmax_rast <- brick("topowx_tmax_2004-2012.grd")
topowx_tmin_rast <- brick("topowx_tmin_2004-2012.grd")

setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/cleaned_dailyair/")
tavg <- read.csv("alldat_daily_tavg_2006_2015.csv")
tavg$DATE <- as.POSIXct(tavg$DATE)
tmax <- read.csv("alldat_daily_tmax_2006_2015.csv")
tmax$DATE <- as.POSIXct(tmax$DATE)
tmin <- read.csv("alldat_daily_tmin_2006_2015.csv")
tmin$DATE <- as.POSIXct(tmin$DATE)
meta <- read.csv("/Users/ian/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/site_metadata.csv")
meta$alt_code <- gsub(pattern="-",replacement=".",x=meta$location,fixed=TRUE)

ibut_date_range <- c(min(tavg$DATE),max(tavg$DATE))

##Drops experiment site outside of study area.
exp_pos <- which(colnames(tavg)=="TAHO.EXP.A1")
tavg <- tavg[,-c(exp_pos)]
tmin <- tmin[,-c(exp_pos)]
tmax <- tmax[,-c(exp_pos)]

#### Data exploration: what does the pattern of spatial dependence look like? ####

##Makes an image plot showing the temperature anomaly at each sensor at each day.
pdf("../results/data_matrix.pdf",width=8,height=10)
ncols <- dim(tavg[,-1])[2]
par(mfrow=c(1,1),mar=c(4,6,1,7),oma=c(1,1,1,1),xpd=TRUE)
image(x=tavg$DATE,y=1:ncols,z=as.matrix(tavg[,-1]),axes=FALSE,
      xlab="",ylab="",col=rainbow(100))
axis(side=2,at=1:ncols,labels=FALSE)
axis.POSIXct(at=seq(as.Date("2007-01-01"),as.Date("2015-01-01"),by="year"),side=1)
image.plot(x=1:dim(tavg)[1],y=1:ncols,z=as.matrix(tavg[,-1]),yaxt="n",
           xlab="",ylab="",legend.only=TRUE)
text(x=as.numeric(min(tavg$DATE))-0.3e7,
     y=1:dim(tavg[,-1])[2],cex=0.6,pos=2,
     labels=colnames(tavg[,-1]))
segments(x0=as.numeric(min(tavg$DATE)),x1=as.numeric(max(tavg$DATE)),
         y0=1:ncols,y1=1:ncols,lty=3,col=rgb(0.2,0.2,0.2,0.2,0.8))
dev.off()

##Plots relationships between predictors.
pdf("../results/covar_pairs.pdf",width=5,height=5)
pairs(~asin(meta$ccov_81m)+meta$elev_9m+meta$relev_729m,
      labels=c("asin(Can. \n Cover)","Elev.\n (m)","Relev.\nElev."),
      cex.labels=1)
dev.off()

#Temporal bounds of the analysis
start_date <- as.POSIXct("2010-7-1")
end_date <- as.POSIXct("2015-7-1")

start_num <- which(tavg$DATE == start_date)
end_num <- which(tavg$DATE == end_date)



##Runs the analysis for tavg, tmax, and tmin.
tavg_diags <- model.temps.lm(tavg,daterange=start_num:end_num)
tmax_diags <- model.temps.lm(tmax,daterange=start_num:end_num)
tmin_diags <- model.temps.lm(tmin,daterange=start_num:end_num)

##Creates plots for each dataset.
pdf("../results/lm_coefficients_diagnostics_tavg.pdf",width=6,height=8)
plot.diags(tavg_diags,maintext="TAVG")
dev.off()
pdf("../results/lm_coefficients_diagnostics_tmax.pdf",width=6,height=8)
plot.diags(tmax_diags,maintext="TMAX")
dev.off()
pdf("../results/lm_coefficients_diagnostics_tmin.pdf",width=6,height=8)
plot.diags(tmin_diags,maintext="TMIN")
dev.off()

#### Plots and measures relationships between observed and regional PRISM and topoWX estimates

##Computes values of disparity, coupling, and buffering in tmax, tmin, and tavg for each site.
indices_tmax <- boot_indices(regional_series=prism$TMAX[976:length(prism$TMAX)],
                             local_frame=tmax,site_names=colnames(tmax[,-1]),nboot=1000)
indices_tmin <- boot_indices(regional_series=prism$TMIN[976:length(prism$TMIN)],
                             local_frame=tmin,site_names=colnames(tmin[,-1]),nboot=1000)
indices_tavg <- boot_indices(regional_series=prism$TEMP[976:length(prism$TEMP)],
                             local_frame=tavg,site_names=colnames(tavg[,-1]),nboot=1000)

## Computes median value for coupling.
med_coup <- median(indices_tmax[,8],na.rm=TRUE)

### Plots each site ranked by disparity,buffering, and coupling.
indices_disprank <- indices_tmax[order(indices_tmax$disp),]
indices_disprank$high <- indices_disprank[,3] > 0
indices_disprank$low <- indices_disprank[,4] < 0
indices_disprank$cat <- as.factor(indices_disprank$low - indices_disprank$high)
indices_sensrank <- indices_tmax[order(indices_tmax$sens),]
indices_sensrank$high <- indices_sensrank[,6] > 1
indices_sensrank$low <- indices_sensrank[,7] < 1 
indices_sensrank$cat <- as.factor(indices_sensrank$low - indices_sensrank$high)
indices_couprank <- indices_tmax[order(indices_tmax$decouple),]
indices_couprank$high <- indices_couprank[,9] > med_coup
indices_couprank$low <- indices_couprank[,10] < med_coup
indices_couprank$cat <- as.factor(indices_couprank$low - indices_couprank$high)

pdf("../results/tmax_disp_sens_coup.pdf",width=8,height=12)
par(mfrow=c(1,3),mar=c(4,1,2,1),xpd=TRUE)
sens_cols <- c("black","black","black")
xmin <- min(indices_disprank[,2],na.rm=TRUE)-5
xmax <- max(indices_disprank[,2],na.rm=TRUE)+3
plot(seq(xmin,xmax,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Disparity (C)",ylab="",main="Disparity",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=-6,x1=2,y0=i,y1=i,lty=3,lwd=0.5)
  points(indices_disprank[i,2],i,pch=20)
  segments(x0=indices_disprank[i,3],x1=indices_disprank[i,4],y0=i,y1=i)
  boxed.labels(labels=indices_disprank[i,1],x=xmin+2,y=i,cex=0.7,border=FALSE,
       col=sens_cols[as.numeric(indices_disprank$cat[i])])
  
}
segments(x0=0,x1=0,y0=-1,y1=i,lty=2)
plot(seq(0.4,1.5,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Sensitivity (C local / C regional)",ylab="",main="Sensitivity",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=0.5,x1=1.5,y0=i,y1=i,lty=3,lwd=0.5)
  points(indices_sensrank[i,5],i,pch=20)
  segments(x0=indices_sensrank[i,6],x1=indices_sensrank[i,7],y0=i,y1=i)
  boxed.labels(labels=indices_sensrank[i,1],x=0.5,y=i,cex=0.7,border=FALSE,
               col=sens_cols[as.numeric(indices_sensrank$cat[i])])
  
}
segments(x0=1,x1=1,y0=-1,y1=i,lty=2)
plot(seq(0,4,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Decoupling (Correlation)",ylab="",main="Decoupling",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=0.88,x1=1,y0=i,y1=i,lty=3,lwd=0.5)
  points(indices_couprank[i,8],i,pch=20)
  segments(x0=indices_couprank[i,9],x1=indices_couprank[i,10],y0=i,y1=i)
  boxed.labels(labels=indices_couprank[i,1],x=0.4,y=i,cex=0.7,
               col=sens_cols[as.numeric(indices_couprank$cat[i])],border=FALSE)
  
}
segments(x0=med_coup,x1=med_coup,y0=-1,y1=i,lty=2)
dev.off()

### Plots each site ranked by disparity,buffering, and coupling.
tmin_indices_disprank <- indices_tmin[order(indices_tmin$disp),]
tmin_indices_disprank$high <- tmin_indices_disprank[,3] > 0
tmin_indices_disprank$low <- tmin_indices_disprank[,4] < 0
tmin_indices_disprank$cat <- as.factor(tmin_indices_disprank$low - tmin_indices_disprank$high)
tmin_indices_sensrank <- indices_tmin[order(indices_tmin$sens),]
tmin_indices_sensrank$high <- tmin_indices_sensrank[,6] > 1
tmin_indices_sensrank$low <- tmin_indices_sensrank[,7] < 1 
tmin_indices_sensrank$cat <- as.factor(tmin_indices_sensrank$low - tmin_indices_sensrank$high)
tmin_indices_couprank <- indices_tmin[order(indices_tmin$decouple),]
tmin_indices_couprank$high <- tmin_indices_couprank[,9] > med_coup
tmin_indices_couprank$low <- tmin_indices_couprank[,10] < med_coup
tmin_indices_couprank$cat <- as.factor(tmin_indices_couprank$low - tmin_indices_couprank$high)

pdf("../results/tmin_disp_sens_coup.pdf",width=8,height=10)
par(mfrow=c(1,3),mar=c(4,1,2,1),xpd=TRUE)
sens_cols <- c("black","black","black")
xmin <- min(tmin_indices_disprank[,2],na.rm=TRUE)-5
xmax <- max(tmin_indices_disprank[,2],na.rm=TRUE)+3
plot(seq(xmin,xmax,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Disparity (C)",ylab="",main="Disparity",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=-6,x1=2,y0=i,y1=i,lty=3,lwd=0.5)
  points(tmin_indices_disprank[i,2],i,pch=20)
  segments(x0=tmin_indices_disprank[i,3],x1=tmin_indices_disprank[i,4],y0=i,y1=i)
  boxed.labels(labels=tmin_indices_disprank[i,1],x=xmin+2,y=i,cex=0.7,border=FALSE,
               col=sens_cols[as.numeric(tmin_indices_disprank$cat[i])])
  
}
segments(x0=0,x1=0,y0=-1,y1=i,lty=2)
plot(seq(0.4,1.5,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Sensitivity (C local / C regional)",ylab="",main="Sensitivity",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=0.5,x1=1.5,y0=i,y1=i,lty=3,lwd=0.5)
  points(tmin_indices_sensrank[i,5],i,pch=20)
  segments(x0=tmin_indices_sensrank[i,6],x1=tmin_indices_sensrank[i,7],y0=i,y1=i)
  boxed.labels(labels=tmin_indices_sensrank[i,1],x=0.5,y=i,cex=0.7,border=FALSE,
               col=sens_cols[as.numeric(tmin_indices_sensrank$cat[i])])
  
}
segments(x0=1,x1=1,y0=-1,y1=i,lty=2)
plot(seq(0,3,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Decoupling (Correlation)",ylab="",main="Decoupling",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=0.88,x1=1,y0=i,y1=i,lty=3,lwd=0.5)
  points(tmin_indices_couprank[i,8],i,pch=20)
  segments(x0=tmin_indices_couprank[i,9],x1=tmin_indices_couprank[i,10],y0=i,y1=i)
  boxed.labels(labels=tmin_indices_couprank[i,1],x=0.3,y=i,cex=0.7,
               col=sens_cols[as.numeric(tmin_indices_couprank$cat[i])],border=FALSE)
  
}
segments(x0=med_coup,x1=med_coup,y0=-1,y1=i,lty=2)
dev.off()

### Plots each site ranked by disparity,buffering, and coupling.
tavg_indices_disprank <- indices_tavg[order(indices_tavg$disp),]
tavg_indices_disprank$high <- tavg_indices_disprank[,3] > 0
tavg_indices_disprank$low <- tavg_indices_disprank[,4] < 0
tavg_indices_disprank$cat <- as.factor(tavg_indices_disprank$low - tavg_indices_disprank$high)
tavg_indices_sensrank <- indices_tavg[order(indices_tavg$sens),]
tavg_indices_sensrank$high <- tavg_indices_sensrank[,6] > 1
tavg_indices_sensrank$low <- tavg_indices_sensrank[,7] < 1 
tavg_indices_sensrank$cat <- as.factor(tavg_indices_sensrank$low - tavg_indices_sensrank$high)
tavg_indices_couprank <- indices_tavg[order(indices_tavg$decouple),]
tavg_indices_couprank$high <- tavg_indices_couprank[,9] > med_coup
tavg_indices_couprank$low <- tavg_indices_couprank[,10] < med_coup
tavg_indices_couprank$cat <- as.factor(tavg_indices_couprank$low - tavg_indices_couprank$high)

pdf("../results/tavg_disp_sens_coup.pdf",width=8,height=10)
par(mfrow=c(1,3),mar=c(4,1,2,1),xpd=TRUE)
sens_cols <- c("black","black","black")
xmin <- min(tavg_indices_disprank[,2],na.rm=TRUE)-5
xmax <- max(tavg_indices_disprank[,2],na.rm=TRUE)+3
plot(seq(xmin,xmax,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Disparity (C)",ylab="",main="Disparity",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=-6,x1=2,y0=i,y1=i,lty=3,lwd=0.5)
  points(tavg_indices_disprank[i,2],i,pch=20)
  segments(x0=tavg_indices_disprank[i,3],x1=tavg_indices_disprank[i,4],y0=i,y1=i)
  boxed.labels(labels=tavg_indices_disprank[i,1],x=xmin+2,y=i,cex=0.7,border=FALSE,
               col=sens_cols[as.numeric(tavg_indices_disprank$cat[i])])
  
}
segments(x0=0,x1=0,y0=-1,y1=i,lty=2)
plot(seq(0.4,1.5,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Sensitivity (C local / C regional)",ylab="",main="Sensitivity",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=0.5,x1=1.5,y0=i,y1=i,lty=3,lwd=0.5)
  points(tavg_indices_sensrank[i,5],i,pch=20)
  segments(x0=tavg_indices_sensrank[i,6],x1=tavg_indices_sensrank[i,7],y0=i,y1=i)
  boxed.labels(labels=tavg_indices_sensrank[i,1],x=0.5,y=i,cex=0.7,border=FALSE,
               col=sens_cols[as.numeric(tavg_indices_sensrank$cat[i])])
  
}
segments(x0=1,x1=1,y0=-1,y1=i,lty=2)
plot(seq(0,3,length.out=10),seq(1,dim(indices_disprank)[1],length.out=10),type="n",
     yaxt="n",xlab="Temp. Decoupling (Correlation)",ylab="",main="Decoupling",frame.plot=FALSE)
for(i in 1:dim(indices_disprank)[1]){
  segments(x0=0.88,x1=1,y0=i,y1=i,lty=3,lwd=0.5)
  points(tavg_indices_couprank[i,8],i,pch=20)
  segments(x0=tavg_indices_couprank[i,9],x1=tavg_indices_couprank[i,10],y0=i,y1=i)
  boxed.labels(labels=tavg_indices_couprank[i,1],x=0.3,y=i,cex=0.7,
               col=sens_cols[as.numeric(tavg_indices_couprank$cat[i])],border=FALSE)
  
}
segments(x0=med_coup,x1=med_coup,y0=-1,y1=i,lty=2)
dev.off()

#### Initial exploration predicting disparity, buffering, and coupling.
meta_indices_tmax <- merge(meta,indices_tmax,by.x="alt_code",by.y="site")
meta_indices_tmin <- merge(meta,indices_tmin,by.x="alt_code",by.y="site")
meta_indices_tavg <- merge(meta,indices_tavg,by.x="alt_code",by.y="site")

###Computes scaled euclidean distance from "best" corner.
best_disp_tavg <- min(scale(meta_indices_tavg$disp),na.rm=TRUE)
best_sens_tavg <- min(scale(meta_indices_tavg$sens),na.rm=TRUE)
best_decoup_tavg <- max(scale(meta_indices_tavg$decouple),na.rm=TRUE)

euc_dist <-  function(disp,sens,decoup){
  sqrt((disp-best_disp_tavg)^2+(sens-best_sens_tavg)^2+(decoup-best_decoup_tavg)^2)
}

dist_tavg <- euc_dist(scale(meta_indices_tavg$disp),
                      scale(meta_indices_tavg$sens),
                      scale(meta_indices_tavg$decouple))
meta_indices_tavg$dist_tavg <- dist_tavg

best_disp_tmin <- min(scale(meta_indices_tmin$disp),na.rm=TRUE)
best_sens_tmin <- min(scale(meta_indices_tmin$sens),na.rm=TRUE)
best_decoup_tmin <- max(scale(meta_indices_tmin$decouple),na.rm=TRUE)

euc_dist <-  function(disp,sens,decoup){
  sqrt((disp-best_disp_tmin)^2+(sens-best_sens_tmin)^2+(decoup-best_decoup_tmin)^2)
}

dist_tmin <- euc_dist(scale(meta_indices_tmin$disp),
                      scale(meta_indices_tmin$sens),
                      scale(meta_indices_tmin$decouple))
meta_indices_tmin$dist_tmin <- dist_tmin

best_disp_tmax <- min(scale(meta_indices_tmax$disp),na.rm=TRUE)
best_sens_tmax <- min(scale(meta_indices_tmax$sens),na.rm=TRUE)
best_decoup_tmax <- max(scale(meta_indices_tmax$decouple),na.rm=TRUE)

euc_dist <-  function(disp,sens,decoup){
  sqrt((disp-best_disp_tmax)^2+(sens-best_sens_tmax)^2+(decoup-best_decoup_tmax)^2)
}

dist_tmax <- euc_dist(scale(meta_indices_tmax$disp),
                      scale(meta_indices_tmax$sens),
                      scale(meta_indices_tmax$decouple))
meta_indices_tmax$dist_tmax <- dist_tmax

####Writes to disk.
write.csv(meta_indices_tavg[!is.na(meta_indices_tavg$decouple),],file="../results/buff_sens_decoup_tavg_covars.csv")
write.csv(meta_indices_tmax[!is.na(meta_indices_tmax$decouple),],file="../results/buff_sens_decoup_tmax_covars.csv")
write.csv(meta_indices_tmin[!is.na(meta_indices_tmin$decouple),],file="../results/buff_sens_decoup_tmin_covars.csv")


## 3-d plot of indices.
require(scatterplot3d)
require(RColorBrewer)

ref_pal <- rev(brewer.pal(n=5,name="Spectral"))
dist_tavg_cl <- cut(meta_indices_tavg$dist_tavg,breaks=5)
dist_tavg_col <- ref_pal[as.numeric(dist_tavg_cl)]

dist_tmin_cl <- cut(meta_indices_tmin$dist_tmin,breaks=5)
dist_tmin_col <- ref_pal[as.numeric(dist_tmin_cl)]

dist_tmax_cl <- cut(meta_indices_tmax$dist_tmax,breaks=5)
dist_tmax_col <- ref_pal[as.numeric(dist_tmax_cl)]


pdf("../results/disp_buff_decoup_3d.pdf",width=4.5,height=12)
par(mfrow=c(3,1))
s3d <- scatterplot3d(x=meta_indices_tavg$disp,y=meta_indices_tavg$sens,z=meta_indices_tavg$decouple,pch=21,
       xlab="Disparity (C)",ylab="Sensitivity",zlab="Decoupling",main="Tavg",
       xlim=c(-8,3),ylim=c(0.7,1.3),zlim=c(0.8,4),color="black",cex.symbols=1.2,bg=dist_tavg_col)
shadow_x <- rep(-8,length(meta_indices_tavg$disp))
s3d$points3d(x=shadow_x,y=meta_indices_tavg$sens,z=meta_indices_tavg$decouple,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)
shadow_y <- rep(1.3,length(meta_indices_tavg$disp))
s3d$points3d(x=meta_indices_tavg$disp,y=shadow_y,z=meta_indices_tavg$decouple,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)
shadow_z <- rep(0.5,length(meta_indices_tavg$disp))
s3d$points3d(x=meta_indices_tavg$disp,y=meta_indices_tavg$sens,z=shadow_z,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)

s3d <- scatterplot3d(x=meta_indices_tmin$disp,y=meta_indices_tmin$sens,z=meta_indices_tmin$decouple,pch=21,
                     xlab="Disparity (C)",ylab="Sensitivity",zlab="Decoupling",main="tmin",
                     xlim=c(-8,3),ylim=c(0.7,1.3),zlim=c(0.8,4),color="black",cex.symbols=1.2,bg=dist_tmin_col)
shadow_x <- rep(-8,length(meta_indices_tmin$disp))
s3d$points3d(x=shadow_x,y=meta_indices_tmin$sens,z=meta_indices_tmin$decouple,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)
shadow_y <- rep(1.3,length(meta_indices_tmin$disp))
s3d$points3d(x=meta_indices_tmin$disp,y=shadow_y,z=meta_indices_tmin$decouple,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)
shadow_z <- rep(0.5,length(meta_indices_tmin$disp))
s3d$points3d(x=meta_indices_tmin$disp,y=meta_indices_tmin$sens,z=shadow_z,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)


s3d <- scatterplot3d(x=meta_indices_tmax$disp,y=meta_indices_tmax$sens,z=meta_indices_tmax$decouple,pch=21,
                     xlab="Disparity (C)",ylab="Sensitivity",zlab="Decoupling",main="tmax",
                     xlim=c(-8,3),ylim=c(0.7,1.3),zlim=c(0.8,4),color="black",cex.symbols=1.2,bg=dist_tmax_col)
shadow_x <- rep(-8,length(meta_indices_tmax$disp))
s3d$points3d(x=shadow_x,y=meta_indices_tmax$sens,z=meta_indices_tmax$decouple,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)
shadow_y <- rep(1.3,length(meta_indices_tmax$disp))
s3d$points3d(x=meta_indices_tmax$disp,y=shadow_y,z=meta_indices_tmax$decouple,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)
shadow_z <- rep(0.5,length(meta_indices_tmax$disp))
s3d$points3d(x=meta_indices_tmax$disp,y=meta_indices_tmax$sens,z=shadow_z,pch=19,
             col=rgb(0.7,0.7,0.7,1),cex.symbols=1.2)

dev.off()

####Plots indices against some obvious covariates.
par(mfrow=c(3,3),xpd=FALSE,mar=c(4,4,2,2),oma=c(1,1,1,1))

##Tmax
plot(disp~elev,pch=20,data=meta_indices_tmax,xlab="Elevation (m)",ylab="Disparity (C)")
segments(x0=meta_indices_tmax$elev,x1=meta_indices_tmax$elev,y0=meta_indices_tmax$disp_lwr95,
         y1=meta_indices_tmax$disp_upr95,lty=1,lwd=0.5)
plot(sens~meta_indices_tmax$srad_sum_9m,pch=20,data=meta_indices_tmax,xlab="Solar Radiation (Wh/m^2/yr)",ylab="Sensitivity (unitless)")
segments(x0=meta_indices_tmax$srad_sum_9m,x1=meta_indices_tmax$srad_sum_9m,y0=meta_indices_tmax$sens_lwr95,
         y1=meta_indices_tmax$sens_upr95,lty=1,lwd=0.5)
plot(decouple~relev_729m,pch=20,data=meta_indices_tmax,xlab="Cold-Air Index",ylab="Decoupling (unitless)")
segments(x0=meta_indices_tmax$relev_729m,x1=meta_indices_tmax$cair_9m,y0=meta_indices_tmax$decouple_lwr95,
         y1=meta_indices_tmax$decouple_upr95,lty=1,lwd=0.5)

##Tmin
plot(disp~elev,pch=20,data=meta_indices_tmin,xlab="Elevation (m)",ylab="Disparity (C)")
segments(x0=meta_indices_tmin$elev,x1=meta_indices_tmin$elev,y0=meta_indices_tmin$disp_lwr95,
         y1=meta_indices_tmin$disp_upr95,lty=1,lwd=0.5)
plot(sens~cvol_81m,pch=20,data=meta_indices_tmin,xlab="Canopy Volume (m^3)",ylab="Sensitivity (unitless)")
segments(x0=meta_indices_tmin$cvol_81m,x1=meta_indices_tmin$cvol_81m,y0=meta_indices_tmin$sens_lwr95,
         y1=meta_indices_tmin$sens_upr95,lty=1,lwd=0.5)
plot(decouple~relev_729m,pch=20,data=meta_indices_tmin,xlab="Rel. Elevation (m)",ylab="Decoupling (unitless)")
segments(x0=meta_indices_tmin$relev_729m,x1=meta_indices_tmin$relev_729m,y0=meta_indices_tmin$decouple_lwr95,
         y1=meta_indices_tmin$decouple_upr95,lty=1,lwd=0.5)

##Tavg
plot(disp~elev,pch=20,data=meta_indices_tavg,xlab="Elevation (m)",ylab="Disparity (C)")
segments(x0=meta_indices_tavg$elev,x1=meta_indices_tavg$elev,y0=meta_indices_tavg$disp_lwr95,
         y1=meta_indices_tavg$disp_upr95,lty=1,lwd=0.5)
plot(sens~cvol_81m,pch=20,data=meta_indices_tavg,xlab="Canopy Volume (m^3)",ylab="Sensitivity (unitless)")
segments(x0=meta_indices_tavg$cvol_81m,x1=meta_indices_tavg$cvol_81m,y0=meta_indices_tavg$sens_lwr95,
         y1=meta_indices_tavg$sens_upr95,lty=1,lwd=0.5)
plot(decouple~relev_729m,pch=20,data=meta_indices_tavg,xlab="Rel. Elevation (m)",ylab="Decoupling (unitless)")
segments(x0=meta_indices_tavg$relev_729m,x1=meta_indices_tavg$relev_729m,y0=meta_indices_tavg$decouple_lwr95,
         y1=meta_indices_tavg$decouple_upr95,lty=1,lwd=0.5)


####Initial models.

disparity_lm_tmax <- lm(disp~elev+cvol_81m,data=meta_indices_tmax)
summary(disparity_lm_tmax)
sens_lm_tmax <- lm(sens~cvol_81m*elev,data=meta_indices_tmax)
summary(sens_lm_tmax)
decoupling_lm_tmax <- lm(decouple~relev_729m+I(relev_729m^2),data=meta_indices_tmax)
summary(decoupling_lm_tmax)

disparity_lm_tmin <- lm(disp~elev+cold_ind_9m+I(logit(ccov_81m+0.001)),data=meta_indices_tmin)
summary(disparity_lm_tmin)
sens_lm_tmin <- lm(sens~elev+I(logit(ccov_81m+0.001))+cold_ind_9m,data=meta_indices_tmin)
summary(sens_lm_tmin)
decoupling_lm_tmin <- lm(decouple~elev+relev_729m,data=meta_indices_tmin)
summary(decoupling_lm_tmin)

disparity_lm_tavg <- lm(disp~elev+cold_ind_9m+I(logit(ccov_81m+0.001)),data=meta_indices_tavg)
summary(disparity_lm_tavg)
sens_lm_tavg <- lm(sens~elev+relev_729m,data=meta_indices_tavg)
summary(sens_lm_tavg)
decoupling_lm_tavg <- lm(decouple~elev,data=meta_indices_tavg)
summary(decoupling_lm_tavg)

#### Fitting a dynamic space-time model to the data.####

##Temporal bounds of analysis
start_date <- as.POSIXct("2013-6-1")
end_date <- as.POSIXct("2015-6-1")

##Subsets data.
tavg_merge <- merge(tavg,prism,all.x=TRUE)
tavg_diff <- tavg_merge[2:dim(tavg)[2]] - tavg_merge$TEMP
tavg_sub <- tavg_diff[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date),-c(11)]
y.tavg <- t(tavg_sub)
colnames(y.tavg) <- paste("y.",1:dim(y.tavg)[2],sep="")
y.prism.tavg <- matrix(tavg_merge$TEMP[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)],
                       nrow=dim(y.tavg)[1],ncol=dim(y.tavg)[2],byrow=TRUE)

tmax_merge <- merge(tmax,prism,all.x=TRUE)
tmax_diff <- tmax_merge[2:dim(tmax)[2]] - tavg_merge$TMAX
tmax_sub <- tmax_diff[which(tmax_merge$DATE == start_date):which(tmax_merge$DATE == end_date),-c(11)]
y.tmax <- t(tmax_sub)
colnames(y.tmax) <- paste("y.",1:dim(y.tmax)[2],sep="")
y.prism.tmax <- matrix(tavg_merge$TMAX[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)],
                       nrow=dim(y.tavg)[1],ncol=dim(y.tavg)[2],byrow=TRUE)

tmin_merge <- merge(tmin,prism,all.x=TRUE)
tmin_diff <- tmin_merge[2:dim(tmin)[2]] - tavg_merge$TMIN
tmin_sub <- tmin_diff[which(tmin_merge$DATE == start_date):which(tmin_merge$DATE == end_date),-c(11)]
y.tmin <- t(tmin_sub)
colnames(y.tmin) <- paste("y.",1:dim(y.tmin)[2],sep="")
y.prism.tmin <- matrix(tavg_merge$TMIN[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)],
                       nrow=dim(y.tavg)[1],ncol=dim(y.tavg)[2],byrow=TRUE)


##Adds Franklin sites to data as missing.
#fsites.y <- matrix(NA,nrow=length(fsites$SITE),ncol=dim(y.tavg)[2])
#rownames(fsites.y) <- fsites$SITE[order(fsites$SITE)]
#colnames(fsites.y) <- paste("y.",1:dim(y.tmin)[2],sep="")
#y.tavg <- rbind(y.tavg,fsites.y)
#y.tmin <- rbind(y.tmin,fsites.y)
#y.tmax <- rbind(y.tmax,fsites.y)

##Dimensions for model fitting and plotting
y.dates <- tavg_merge$DATE[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)]
y.prcp <- tavg_merge$PRCP[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)]
tavg.prsm <- tavg_merge$TEMP[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)]
tmin.prsm <- tmin_merge$TMIN[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)]
tmax.prsm <- tmax_merge$TMAX[which(tavg_merge$DATE == start_date):which(tavg_merge$DATE == end_date)]
N.t<- dim(tavg_sub)[1]
n <- nrow(y.tavg)

##Drops sites with no covariate data and assembles covariates.
cols <- data.frame(alt_code=colnames(tavg_sub))
meta_order <-merge(cols,meta,by.x="alt_code",by.y="alt_code",sort=FALSE,all.x=TRUE)

##Prepares covariates.
meta_covs <- data.frame(site=meta_order$alt_code,
                        elev=meta_order$elev,
                        relev=meta_order$relev_729m,
                        cvol=meta_order$cvol_81m,
                        ccov=meta_order$ccov_81m,
                        cair=meta_order$cold_ind_9m,
                        srad=meta_order$srad_sum_9m,
                        dry=meta_order$dry_ind_9m,
                        utmx=meta_order$utmx,
                        utmy=meta_order$utmy)
# f_covs <- data.frame(site=fsites$SITE,
#                      elev=fsites$MORA_elev_3m,
#                      relev=(fsites$MORA_elev_3m - fsites$MORA_elev__focal_729m),
#                      cvol=fsites$MORA_can_vol_81m,
#                      ccov=fsites$MORA_can_pct_focal81m,
#                      cair=fsites$MORA_coldair_index,
#                      utmx=fsites$utmx,
#                      utmy=fsites$utmy)
#meta_covs <- rbind(meta_covs,f_covs)


elev <- scale(meta_covs$elev/1000)
relev <- scale(meta_covs$relev)
cvol <- scale(sqrt(meta_covs$cvol))
cvol_resid <- scale(lm(cvol~elev)$residuals)
ccov1 <- meta_covs$ccov
ccov1[ccov1<0.01] <- 0.01
ccov1[ccov1>0.99] <- 0.99
ccov <- scale(logit(ccov1))
ccov_resid <- scale(lm(ccov~elev)$residuals)
srad <- scale(meta_covs$srad)
cair <- scale(meta_covs$cair)
cair_resid <- scale(lm(cair~elev)$residuals)
srad <- scale(meta_covs$srad)
srad_resid <- scale(lm(srad~elev)$residuals)
dry <- scale(meta_covs$dry)
dry_resid <- scale(lm(dry~elev)$residuals)
utmx <- scale(meta_covs$utmx)
utmy <- scale(meta_covs$utmy)

##Gets coordinates in km.
coords1 <- as.matrix(meta_covs[,c("utmx","utmy")]/1000)
coords <- jitter(coords1,amount=0.001)
max.d <- max(iDist(coords))

##make symbolic model formula statement for each day
mods_tavg <- lapply(paste(colnames(y.tavg),'~','elev*cair_resid + cvol_resid + utmx + utmy',sep=''), as.formula)
mods_tmin <- lapply(paste(colnames(y.tmin),'~','elev*cair_resid + cvol_resid + utmx + utmy',sep=''), as.formula)
mods_tmax <- lapply(paste(colnames(y.tmax),'~','elev*cair_resid + cvol_resid + utmx + utmy',sep=''), as.formula)

p <- 7 # num. linear terms

##add some missing observations to illustrate prediction
# miss <- sample(1:N.t, 15)
# holdout.station.id<- 20
# y.t.holdout <- y.t[holdout.station.id, miss]
# y.t[holdout.station.id, miss] <- NA

##set starting and priors

starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
                 "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                 "sigma.eta"=diag(rep(0.01, p)))

tuning <- list("phi"=rep(5, N.t)) 

priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
               "phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
               "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
               "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
               "sigma.eta.IW"=list(2, diag(0.001,p)))

n.samples <- 2000

m.1 <- spDynLM(mods_tavg, data=data.frame(cbind(y.tavg,elev,cair,cvol_resid,utmx,utmy)), 
                                                   coords=coords,starting=starting, 
                                                   tuning=tuning, priors=priors, 
                                                   get.fitted =TRUE,cov.model="spherical", 
                                                   n.samples=n.samples, n.report=10)
save(m.1,file="/volumes/ib_working/tavg_ssmod_output.Rdata",compress=TRUE)
remove(m.1)

m.2 <- spDynLM(mods_tmax, data=data.frame(cbind(y.tmax,elev,cair,cvol_resid,utmx,utmy)), 
                                                   coords=coords,starting=starting, 
                                                   tuning=tuning, priors=priors, 
                                                   get.fitted =TRUE,cov.model="spherical", 
                                                   n.samples=n.samples, n.report=10)
save(m.2,file="/volumes/ib_working/tmax_ssmod_output.Rdata",compress=TRUE)
remove(m.2)

m.3 <- spDynLM(mods_tmin, data=data.frame(cbind(y.tavg,elev,cair,cvol_resid,utmx,utmy)), 
                                            coords=coords,starting=starting, 
                                            tuning=tuning, priors=priors, 
                                            get.fitted =TRUE,cov.model="spherical", 
                                            n.samples=n.samples, n.report=10)
save(m.3,file="/volumes/ib_working/tmin_ssmod_output.Rdata",compress=TRUE)
remove(m.3)

load("/volumes/ib_working/tavg_ssmod_output.Rdata")
load("/volumes/ib_working/tmax_ssmod_output.Rdata")
load("/volumes/ib_working/tmin_ssmod_output.Rdata")

burn.in <- floor(0.75*n.samples)

quant <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}

## Binds all coefficients into a data frame
beta <- apply(m.1$p.beta.samples[burn.in:n.samples,], 2, quant)
beta.0 <- beta[,grep("Intercept", colnames(beta))]
beta.1 <- beta[,grep("^elev\\.", colnames(beta))]
beta.2 <- beta[,grep("^cair_resid\\.", colnames(beta))]
beta.3 <- beta[,grep("elev:cair_resid", colnames(beta))]
beta.4 <- beta[,grep("cvol", colnames(beta))]
beta.5 <- beta[,grep("^utmx\\.", colnames(beta))]
beta.6 <- beta[,grep("^utmy\\.", colnames(beta))]
tavg.betas <- data.frame(y.dates,beta.0[1,],beta.1[1,],beta.2[1,],beta.3[1,],beta.4[1,],beta.5[1,],beta.6[1,])
tavg.betas.lwr <- data.frame(y.dates,beta.0[2,],beta.1[2,],beta.2[2,],beta.3[2,],beta.4[2,],beta.5[2,],beta.6[2,])
tavg.betas.upr <- data.frame(y.dates,beta.0[3,],beta.1[3,],beta.2[3,],beta.3[3,],beta.4[3,],beta.5[3,],beta.6[3,])
tavg.betas$precip <- y.prcp
tavg.betas$prsm <- tavg.prsm
colnames(tavg.betas) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy","prsm_prcp","prsm_temp")
rownames(tavg.betas) <- y.dates
colnames(tavg.betas.lwr) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy")
rownames(tavg.betas.lwr) <- y.dates
colnames(tavg.betas.upr) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy")
rownames(tavg.betas.upr) <- y.dates

beta.tmax <- apply(m.2$p.beta.samples[burn.in:n.samples,], 2, quant)
beta.tmax.0 <- beta.tmax[,grep("Intercept", colnames(beta.tmax))]
beta.tmax.1 <- beta.tmax[,grep("^elev\\.", colnames(beta.tmax))]
beta.tmax.2 <- beta.tmax[,grep("^cair_resid\\.", colnames(beta.tmax))]
beta.tmax.3 <- beta.tmax[,grep("elev:cair_resid", colnames(beta.tmax))]
beta.tmax.4 <- beta.tmax[,grep("cvol", colnames(beta.tmax))]
beta.tmax.5 <- beta.tmax[,grep("^utmx\\.", colnames(beta.tmax))]
beta.tmax.6 <- beta.tmax[,grep("^utmy\\.", colnames(beta.tmax))]
tmax.betas <- data.frame(y.dates,beta.tmax.0[1,],beta.tmax.1[1,],beta.tmax.2[1,],
                         beta.tmax.3[1,],beta.tmax.4[1,],beta.tmax.5[1,],beta.tmax.6[1,])
tmax.betas.lwr <- data.frame(y.dates,beta.tmax.0[2,],beta.tmax.1[2,],beta.tmax.2[2,],
                             beta.tmax.3[2,],beta.tmax.4[2,],beta.tmax.5[2,],beta.tmax.6[2,])
tmax.betas.upr <- data.frame(y.dates,beta.tmax.0[3,],beta.tmax.1[3,],beta.tmax.2[3,],beta.tmax.3[3,],
                             beta.tmax.4[3,],beta.tmax.5[3,],beta.tmax.6[3,])
tmax.betas$precip <- y.prcp
tmax.betas$prsm <- tmax.prsm
colnames(tmax.betas) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy","prsm_prcp","prsm_temp")
rownames(tmax.betas) <- y.dates
colnames(tmax.betas.lwr) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy")
rownames(tmax.betas.lwr) <- y.dates
colnames(tmax.betas.upr) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy")
rownames(tmax.betas.upr) <- y.dates

beta.tmin <- apply(m.3$p.beta.samples[burn.in:n.samples,], 2, quant)
beta.tmin.0 <- beta.tmin[,grep("Intercept", colnames(beta.tmin))]
beta.tmin.1 <- beta.tmin[,grep("^elev\\.", colnames(beta.tmin))]
beta.tmin.2 <- beta.tmin[,grep("^cair_resid\\.", colnames(beta.tmin))]
beta.tmin.3 <- beta.tmin[,grep("elev:cair_resid", colnames(beta.tmin))]
beta.tmin.4 <- beta.tmin[,grep("cvol", colnames(beta.tmin))]
beta.tmin.5 <- beta.tmin[,grep("^utmx\\.", colnames(beta.tmin))]
beta.tmin.6 <- beta.tmin[,grep("^utmy\\.", colnames(beta.tmin))]
tmin.betas <- data.frame(y.dates,beta.tmin.0[1,],beta.tmin.1[1,],beta.tmin.2[1,],
                         beta.tmin.3[1,],beta.tmin.4[1,],beta.tmin.5[1,],beta.tmin.6[1,])
tmin.betas.lwr <- data.frame(y.dates,beta.tmin.0[2,],beta.tmin.1[2,],beta.tmin.2[2,],
                             beta.tmin.3[2,],beta.tmin.4[2,],beta.tmin.5[2,],beta.tmin.6[2,])
tmin.betas.upr <- data.frame(y.dates,beta.tmin.0[3,],beta.tmin.1[3,],beta.tmin.2[3,],beta.tmin.3[3,],
                             beta.tmin.4[3,],beta.tmin.5[3,],beta.tmin.6[3,])
tmin.betas$precip <- y.prcp
tmin.betas$prsm <- tmin.prsm
colnames(tmin.betas) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy","prsm_prcp","prsm_temp")
rownames(tmin.betas) <- y.dates
colnames(tmin.betas.lwr) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy")
rownames(tmin.betas.lwr) <- y.dates
colnames(tmin.betas.upr) <-  c("date","int","elev","cair","elev_cair","cvol","utmx","utmy")
rownames(tmin.betas.upr) <- y.dates

##Plots time-series of coefficients
plot.beta.ts(y.dates,beta.0,beta.1,beta.2,beta.3,beta.4)
plot.beta.ts(y.dates,beta.tmax.0,beta.tmax.1,beta.tmax.2,beta.tmax.3,beta.tmax.4)
plot.beta.ts(y.dates,beta.tmin.0,beta.tmin.1,beta.tmin.2,beta.tmin.3,beta.tmax.4)

##Function to determine which days have coefficients credibly greater than zero.
cred_cat_fun <- function(lwr_frame,upr_frame) {
  above_yn <- lwr_frame[,-1] > 0
  below_yn <- upr_frame[,-1] < 0
  out <- above_yn
  out[,] <- "zero"
  out[above_yn] <- "above"
  out[below_yn] <- "below"
  colnames(out) <- paste("cred_",colnames(out),sep="")
  return(out)
}

tavg.betas.cred <- cred_cat_fun(tavg.betas.lwr,tavg.betas.upr)
tmin.betas.cred <- cred_cat_fun(tmin.betas.lwr,tmin.betas.upr)
tmax.betas.cred <- cred_cat_fun(tmax.betas.lwr,tmax.betas.upr)

tavg.betas.all <- cbind(tavg.betas,tavg.betas.cred)
tmin.betas.all <- cbind(tmin.betas,tmin.betas.cred)
tmax.betas.all <- cbind(tmax.betas,tmax.betas.cred)

## Extracts posterior samples of coefficients by season.
t.month <- as.numeric(format(y.dates,format="%m"))
t.summer <- t.month %in% c(6,7,8)
t.fall <- t.month %in% c(9,10,11)
t.winter <- t.month %in% c(12,1,2)
t.spring <- t.month %in% c(3,4,5)
t.gs <- t.month %in% c(5,6,7,8,9)

tavg.betas.gs <- tavg.betas.all[t.gs,]
tavg.betas.gs$meas <- "tavg"
tavg.betas.gs$seas <- "growing"
tmax.betas.winter <- tmax.betas.all[t.winter,]
tmax.betas.winter$meas <- "tmax"
tmax.betas.winter$seas <- "winter"
tmin.betas.winter <- tmin.betas.all[t.winter,]
tmin.betas.winter$meas <- "tmin"
tmin.betas.winter$seas <- "winter"
tmax.betas.summer <- tmax.betas.all[t.summer,]
tmax.betas.summer$meas <- "tmax"
tmax.betas.summer$seas <- "summer"
tmin.betas.summer <- tmin.betas.all[t.summer,]
tmin.betas.summer$meas <- "tmin"
tmin.betas.summer$seas <- "summer"

## Binds all coefficients into a data frame for plotting.
require(reshape2)
require(plyr)
betas.all <- rbind(tavg.betas.gs,tmax.betas.winter,tmin.betas.winter,tmax.betas.summer,tmin.betas.summer)
betas.all$seas_meas <- as.factor(paste(betas.all$seas,betas.all$meas))
betas.all$prsm_prcp <- as.factor(betas.all$prsm_prcp > 0.01)
betas.long <- melt(betas.all,measure.vars=c("int","elev","cair","elev_cair","cvol"))
betas.long$variable <- factor(betas.long$variable,labels=c("Intercept","Elevation",
                                                           "Cold-Air","Elev.:Cold-Air","Canopy Volume"))


## Computes the percentage of days with credible coefficients for each season.
betas.cred.pct <- ddply(betas.long,.(seas_meas,prsm_prcp),summarize,
                        int.perc_above = sum(cred_int=="above") / length(cred_int),
                        int.perc_zero = sum(cred_int=="zero") / length(cred_int),
                        int.perc_below = sum(cred_int=="below") / length(cred_int),
                        elev.perc_above = sum(cred_elev=="above") / length(cred_int),
                        elev.perc_zero = sum(cred_elev=="zero") / length(cred_int),
                        elev.perc_below = sum(cred_elev=="below") / length(cred_int),
                        cair.perc_above = sum(cred_cair=="above") / length(cred_int),
                        cair.perc_zero = sum(cred_cair=="zero") / length(cred_int),
                        cair.perc_below = sum(cred_cair=="below") / length(cred_int),
                        elev_cair.perc_above = sum(cred_elev_cair=="above") / length(cred_int),
                        elev_cair.perc_zero = sum(cred_elev_cair=="zero") / length(cred_int),
                        elev_cair.perc_below = sum(cred_elev_cair=="below") / length(cred_int))
betas.cred.pct.long <- melt(betas.cred.pct,id.vars=c("seas_meas","prsm_prcp"))
betas.cred.pct.long$parameter <- matrix(unlist(strsplit(as.character(betas.cred.pct.long$variable),split="\\.")),
                                        nrow=length(betas.cred.pct.long$variable),byrow=TRUE)[,1]
betas.cred.pct.long$stat <- matrix(unlist(strsplit(as.character(betas.cred.pct.long$variable),split="\\.")),
                                   nrow=length(betas.cred.pct.long$variable),byrow=TRUE)[,2]

betas.cred.pct.long$stat <- factor(betas.cred.pct.long$stat,levels=c("perc_above","perc_zero","perc_below"),
                                   labels=c("> 0","- 0 -","< 0"))
betas.cred.pct.long$parameter <- factor(betas.cred.pct.long$parameter,levels=c("int","elev","cair","elev_cair","cvol"),
                                        labels=c("Intercept","Elevation","Cold-Air","Elev.:Cold-Air","Canopy Volume"))
betas.cred.pct.long$seas_meas <- factor(betas.cred.pct.long$seas_meas,
                                        labels=c("G.S. Tavg","Summer Tmax","Summer Tmin","Winter Tmax,","Winter Tmin"))

##Computes median parameters for each season and measurement.
## Computes the percentage of days with credible coefficients for each season.
betas.cred.med <- ddply(betas.long,.(variable,seas_meas,prsm_prcp),summarize,
                        med =  quantile(value,probs=c(0.5)),
                        lwr_quart = quantile(value,probs=c(0.25)),
                        upr_quart = quantile(value,probs=c(0.75)))
prsm.cred.med <- ddply(betas.long,.(seas_meas,prsm_prcp),summarize,
                        temp_med =  quantile(prsm_temp,probs=c(0.5)),
                        temp_lwr_quart = quantile(prsm_temp,probs=c(0.25)),
                        temp_upr_quart = quantile(prsm_temp,probs=c(0.75)))
betas.med.wide <- dcast(betas.cred.med,formula=seas_meas+prsm_prcp~variable,value.var = "med")
betas.med.wide$prsm_temp_med <- prsm.cred.med$temp_med
betas.med.wide$seas_meas_prcp <- paste(betas.med.wide$seas_meas,betas.med.wide$prsm_prcp,sep=" ")

##Predictions from median coefficient estimates.
newelev <- seq(range(elev)[1],range(elev)[2],length.out=1000)
newelev_unscaled <- (newelev * attr(elev,"scaled:scale") + attr(elev,"scaled:center"))*1000

newcair1 <- rep(quantile(cair,probs=0.25),length(newelev))
newcair2 <- rep(quantile(cair,probs=0.75),length(newelev))

reg_fun_ridge <- function(x){x[1] + x[2]*newelev + x[3]*newcair1 + x[4]*newelev*newcair1}
reg_fun_cove <- function(x){x[1] + x[2]*newelev + x[3]*newcair2 + x[4]*newelev*newcair2}

pred_ridge <- apply(betas.med.wide[,3:6],MARGIN=1,FUN=reg_fun_ridge)
colnames(pred_ridge) <- paste("ridge",betas.med.wide$seas_meas,betas.med.wide$prsm_prcp,sep=" ")
rownames(pred_ridge) <- 1:dim(pred_ridge)[1]
pred_cove <- apply(betas.med.wide[,3:6],MARGIN=1,FUN=reg_fun_cove)
colnames(pred_cove) <- paste("cove",betas.med.wide$seas_meas,betas.med.wide$prsm_prcp,sep=" ")
rownames(pred_cove) <- 1:dim(pred_cove)[1]
preds <- data.frame(cbind(pred_ridge,pred_cove))
preds.long <- melt(preds)
vars <- strsplit(as.character(preds.long$variable),split="\\.")
vars_df <- data.frame(do.call(rbind, vars))
colnames(vars_df) <- c("Topo","Season","Meas","Precip")
preds_vars <- data.frame(vars_df,preds.long,elev=rep(newelev_unscaled,
                                                     dim(preds.long)[1]/length(newelev_unscaled)))
preds_vars$seas_meas <- paste(preds_vars$Season,preds_vars$Meas,sep=" ")
  

##Stacked bar chart for credibiltiy intervals.
require(ggplot2)
a <- ggplot(data=subset(betas.cred.pct.long,parameter!="Intercept"))+
     geom_bar(aes(x=parameter,y=value,alpha=stat,fill=parameter,
                  order= -as.numeric(stat)),stat="identity",width=0.4)+
     scale_y_continuous(breaks=c(0,0.5,1))+
     scale_alpha_discrete(breaks=c(0.2,0.6,0.9))+
     xlab("")+
     ylab("Proportion Credible")+
     facet_grid(prsm_prcp~seas_meas)+
     theme_bw()+
     theme( axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.margin = unit(c(2,2,0,2),"mm"))

##Violin plot of seasonal estimates
b <- ggplot(data=subset(betas.long,variable!="Intercept"))+
     geom_violin(aes(x=variable,y=value,fill=variable),linetype="blank")+
     geom_hline(aes(y=0),linetype=2)+
     facet_grid(prsm_prcp~seas_meas)+
     xlab("")+
     ylab("Coefficient Value")+
     theme_bw()+
     theme( axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.margin = unit(c(2,2,0,2),"mm"))

c <- ggplot(preds_vars)+
        geom_line(aes(x=elev,y=value,linetype=Topo,color=Topo),
                  lwd=1.5)+
        scale_color_grey(start=0.2,end=0.5)+
        facet_grid(Precip~seas_meas)+
        xlab("Elevation (m)")+
        ylab("Temperature anomaly (C)")+
        theme_bw()+
        theme(panel.grid.major.x = element_blank(),
              plot.margin = unit(c(2,10,0,2),"mm"))

plotlist <- list(a,b,c)
pdf("../results/ssmod_coefs_violin_2013_2015.pdf",width=8,height=8)
multiplot(plotlist=plotlist)
dev.off()

theta <- apply(m.1$p.theta.samples[burn.in:n.samples,], 2, quant)
sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
tau.sq <- theta[,grep("tau.sq", colnames(theta))]
phi <- theta[,grep("phi", colnames(theta))]

par(mfrow=c(3,1))
plot(1:N.t, sigma.sq[1,], pch=19, cex=0.5, xlab="months", ylab="sigma.sq", ylim=range(sigma.sq))
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[3,], length=0.02, angle=90)
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[2,], length=0.02, angle=90)

plot(1:N.t, tau.sq[1,], pch=19, cex=0.5, xlab="months", ylab="tau.sq", ylim=range(tau.sq))
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[3,], length=0.02, angle=90)
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[2,], length=0.02, angle=90)

plot(1:N.t, 3/phi[1,], pch=19, cex=0.5, xlab="months", ylab="eff. range (km)", ylim=range(3/phi))
arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[3,], length=0.02, angle=90)
arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[2,], length=0.02, angle=90)

y.hat.tavg <- apply(m.1$p.y.samples[,burn.in:n.samples], 1, quant)
y.hat.tavg.med <- matrix(y.hat.tavg[1,], ncol=N.t) + y.prism.tavg
y.hat.tavg.up <- matrix(y.hat.tavg[3,], ncol=N.t) + y.prism.tavg
y.hat.tavg.low <- matrix(y.hat.tavg[2,], ncol=N.t) + y.prism.tavg

y.hat.tmax <- apply(m.2$p.y.samples[,burn.in:n.samples], 1, quant)
y.hat.tmax.med <- matrix(y.hat.tmax[1,], ncol=N.t) + y.prism.tmax
y.hat.tmax.up <- matrix(y.hat.tmax[3,], ncol=N.t) + y.prism.tmax
y.hat.tmax.low <- matrix(y.hat.tmax[2,], ncol=N.t) + y.prism.tmax

y.hat.tmin <- apply(m.3$p.y.samples[,burn.in:n.samples], 1, quant)
y.hat.tmin.med <- matrix(y.hat.tmin[1,], ncol=N.t) + y.prism.tmin
y.hat.tmin.up <- matrix(y.hat.tmin[3,], ncol=N.t) + y.prism.tmin
y.hat.tmin.low <- matrix(y.hat.tmin[2,], ncol=N.t) + y.prism.tmin

pdf("../results/ssmod_2014_series.pdf",width=7,height=4.5)
plots <- c("TO04.STR.A1","TO04.A2","PARA.STR.A1","PARA.A2")
par(mfrow=c(length(plots),1),mar=c(0,2,0,0),oma=c(6,3,3,1))
for(i in 1:length(plots)){
  plotnum <- which(rownames(y.tavg) == plots[i])
  plot(y.dates,(y.tavg[plotnum,]+y.prism.tavg[1,]), pch=19, cex=0.2, xlab="",
       ylab="",ylim=c(-24,34),xaxt="n")
  if(i==1){
    legend(x=as.numeric(y.dates[100]),y=60,title="Modeled",xpd=NA,legend=c("Tmin","Tavg","Tmax"),
           pt.cex=0.1,ncol=3,bty="n",border=rgb(1,1,1,0),
           fill=c(rgb(0.25,0.25,0.9,0.25),rgb(0.25,0.25,0.25,0.25),rgb(0.9,0.25,0.25,0.25)))
  }
  if(i==1){
    legend(x=as.numeric(y.dates[400]),y=60,title="Measured",xpd=NA,legend=c("Tmin","Tavg","Tmax"),pch=c(1,1,1),
           pt.cex=0.4,ncol=3,bty="n",col=c("blue","black","red"))
  }
  polygon(x=c(y.dates,rev(y.dates)),y=c(y.hat.tavg.low[plotnum,],rev(y.hat.tavg.up[plotnum,])),
              col=rgb(0.25,0.25,0.25,0.25),border = FALSE)
  polygon(x=c(y.dates,rev(y.dates)),y=c(y.hat.tmax.low[plotnum,],rev(y.hat.tmax.up[plotnum,])),
          col=rgb(0.9,0.25,0.25,0.25),border = FALSE)
  polygon(x=c(y.dates,rev(y.dates)),y=c(y.hat.tmin.low[plotnum,],rev(y.hat.tmin.up[plotnum,])),
          col=rgb(0.25,0.25,0.9,0.25),border = FALSE)
  points(y.dates,(y.tmax[plotnum,]+y.prism.tmax[1,]), pch=19, cex=0.3, xlab="",
         ylab="",ylim=c(-8,28),xaxt="n",col="red")
  points(y.dates,(y.tmin[plotnum,]+y.prism.tmin[1,]), pch=19, cex=0.3, xlab="",
         ylab="",ylim=c(-8,28),xaxt="n",col="blue")
  axis.POSIXct(side=1,x = y.dates,labels = FALSE)
  abline(h=0,lty=2,col="blue")
  mtext(text=plots[i],side=3,outer=FALSE,padj=2,adj=1)
}
axis.POSIXct(side=1, x=y.dates,labels=TRUE)
mtext(text="Daily Air Temperature (C)",side = 2,outer = TRUE,padj=-1)
mtext(text="Year",side = 1,outer = TRUE,padj=4)
dev.off()

#### Prediction: Using the posterior estimates to predict data at new locations.####
