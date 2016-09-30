##Script to compute figures and summaries of MORAclim data.

##Sets up workspace.
library(raster)
library(maptools)
library(ggplot2)

##Loads data.
ann_tavg <- brick("/Volumes/ib_working/MORAclim_summaries/ann_tavg_sum.tif")
ann_tmin <- brick("/Volumes/ib_working/MORAclim_summaries/ann_tmin_sum.tif")
ann_tmax <- brick("/Volumes/ib_working/MORAclim_summaries/ann_tmax_sum.tif")

sum_tmax <- brick("/Volumes/ib_working/MORAclim_summaries/jja_tmax_sum.tif")
win_tmin <- brick("/Volumes/ib_working/MORAclim_summaries/djf_tmin_sum.tif")

spring_tmax <- brick("/Volumes/ib_working/MORAclim_summaries/mam_tmax_sum.tif")
spring_tmin <- brick("/Volumes/ib_working/MORAclim_summaries/mam_tmin_sum.tif")

sfgdd <- raster("/Volumes/ib_working/MORAclim_summaries/sfgdd_avg_2010_2015.tif")
sffdd <- raster("/Volumes/ib_working/MORAclim_summaries/sffdd_avg_2010_2015.tif")

mora_bound <- readShapePoly("/Volumes/ib_working/GIS/MORA_boundary.shp")

##Resamples predictors to 90m resolution.
preds <- brick("/Volumes/ib_working/MORAclim_preds/moraclim_preds_2016.grd")
preds_90m <- brick("~/GIS/mcpreds_2016_90m.grd")
# preds_90m <- aggregate(preds,fact=30,expand=TRUE,na.rm=TRUE,
#                        filename="~/GIS/mcpreds_2016_90m.grd")

##Creates a map with MAT deviation from elevation trend.
mat_data <- data.frame(tavg=ann_tavg[[1]][],elev=preds_90m[[5]][])
mat_lm <- lm(tavg~elev,data=mat_data)
mat_pred <- preds_90m[[1]]
mat_pred[] <- predict(mat_lm,newdata=mat_data)
mat_resid <- ann_tavg[[1]] - mat_pred
writeRaster(mat_resid,filename="~/GIS/ann_tavg_resid.tif",overwrite=TRUE)

##Creates a map with summer tmax deviation from elevation trend.
sum_tmax_data <- data.frame(tmax=sum_tmax[[1]][],elev=preds_90m[[5]][])
sum_tmax_lm <- lm(tmax~elev,data=sum_tmax_data)
sum_tmax_pred <- preds_90m[[1]]
sum_tmax_pred[] <- predict(sum_tmax_lm,newdata=sum_tmax_data)
sum_tmax_resid <- sum_tmax[[1]] - sum_tmax_pred
writeRaster(sum_tmax_resid,filename="~/GIS/sum_tmax_resid.tif",overwrite=TRUE)


##Creates a map with extreme summer tmax deviation from elevation trend.
ext_tmax_data <- data.frame(tmax=sum_tmax[[5]][],elev=preds_90m[[5]][])
ext_tmax_lm <- lm(tmax~elev,data=ext_tmax_data)
ext_tmax_pred <- preds_90m[[1]]
ext_tmax_pred[] <- predict(ext_tmax_lm,newdata=ext_tmax_data)
ext_tmax_resid <- sum_tmax[[5]] - ext_tmax_pred
writeRaster(ext_tmax_resid,filename="~/GIS/ext_tmax_resid.tif",overwrite=TRUE)


##Creates a map with winter tmin deviation from elevation trend.
win_tmin_data <- data.frame(tmin=win_tmin[[1]][],elev=preds_90m[[5]][])
win_tmin_lm <- lm(tmin~elev,data=win_tmin_data)
win_tmin_pred <- preds_90m[[1]]
win_tmin_pred[] <- predict(win_tmin_lm,newdata=win_tmin_data)
win_tmin_resid <- win_tmin[[1]] - win_tmin_pred
writeRaster(win_tmin_resid,filename="~/GIS/win_tmin_resid.tif",overwrite=TRUE)


##Creates a map with extreme winter tmin deviation from elevation trend.
ext_tmin_data <- data.frame(tmin=win_tmin[[2]][],elev=preds_90m[[5]][])
ext_tmin_lm <- lm(tmin~elev,data=ext_tmin_data)
ext_tmin_pred <- preds_90m[[1]]
ext_tmin_pred[] <- predict(ext_tmin_lm,newdata=ext_tmin_data)
ext_tmin_resid <- win_tmin[[2]] - ext_tmin_pred
writeRaster(ext_tmin_resid,filename="~/GIS/ext_tmin_resid.tif",overwrite=TRUE)

##Creates a map with snow-free growing-degree-days deviation from elevation trend.
sfgdd_data <- data.frame(sfgdd=sfgdd[],elev=preds_90m[[5]][])
sfgdd_lm <- lm(sfgdd~elev,data=sfgdd_data)
sfgdd_pred <- preds_90m[[1]]
sfgdd_pred[] <- predict(sfgdd_lm,newdata=sfgdd_data)
sfgdd_resid <- sfgdd - sfgdd_pred
writeRaster(sfgdd_resid,filename="~/GIS/sfgdd_resid.tif",overwrite=TRUE)


##Creates a map with extreme temperature ranges
ann_range <- sum_tmax[[5]] - win_tmin[[2]]

##Creates plots

setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/")
pdf("MAT_map.pdf",width=5,height=6)
par(mar=c(0.5,0.5,0.5,0.5))
plot(ann_tavg[[1]],col=rev(rainbow(n=30)),box=FALSE,axes=FALSE)
plot(mora_bound,lwd=2,add=TRUE)
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
title("Mean Annual. Temp.",line=-4)
dev.off()

setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/")
pdf("resid_MAT_map.pdf",width=5,height=6)
par(mar=c(0.5,0.5,0.5,0.5))
resid_pal <- colorRampPalette(colors=c("darkblue","white","red"))
resid_seq <- c(-6,-4,-2,-1,-0.5,0,0.5,1,2,4,6)
plot(mat_resid,col=resid_pal(n=length(resid_seq)),breaks=resid_seq,box=FALSE,axes=FALSE)
plot(mora_bound,lwd=2,add=TRUE)
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
title("Mean Annual. Temp.",line=-4)
dev.off()

pdf("resid_sum_tmax_map.pdf",width=5,height=6)
par(mar=c(0.5,0.5,0.5,0.5))
plot(sum_tmax_resid,col=resid_pal(n=length(resid_seq)),breaks=resid_seq,box=FALSE,axes=FALSE)
plot(mora_bound,lwd=2,add=TRUE)
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
title("Summer Tmax",line=-4)
dev.off()

pdf("resid_win_tmin_map.pdf",width=5,height=6)
par(mar=c(0.5,0.5,0.5,0.5))
plot(win_tmin_resid,col=resid_pal(n=length(resid_seq)),breaks=resid_seq,box=FALSE,axes=FALSE)
plot(mora_bound,lwd=2,add=TRUE)
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
title("Winter Tmin",line=-4)
dev.off()

pdf("resid_ext_tmax_map.pdf",width=5,height=6)
par(mar=c(0.5,0.5,0.5,0.5))
plot(ext_tmax_resid,col=resid_pal(n=length(resid_seq)),breaks=resid_seq,box=FALSE,axes=FALSE)
plot(mora_bound,lwd=2,add=TRUE)
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
title("Extreme Tmax",line=-4)
dev.off()

pdf("resid_ext_tmin_map.pdf",width=5,height=6)
par(mar=c(0.5,0.5,0.5,0.5))
plot(ext_tmin_resid,col=resid_pal(n=length(resid_seq)),breaks=resid_seq,box=FALSE,axes=FALSE)
plot(mora_bound,lwd=2,add=TRUE)
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
title("Extreme Tmin",line=-4)
dev.off()

##Plots annual temperature ranges.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/results/")
pdf("ext_range_map.pdf",width=5,height=6)
par(mar=c(0.5,0.5,0.5,0.5))
ext_pal <- colorRampPalette(colors=c("white","yellow","orange","red"))
ext_seq <- seq(28,48,by=2)
plot(ann_range,col=ext_pal(n=length(ext_seq)),breaks=ext_seq,box=FALSE,axes=FALSE)
plot(mora_bound,lwd=2,add=TRUE)
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
title("Extreme Temp. Range",line=-4)
dev.off()

##Clusters climate residuals using kmeans.
clust_data <- data.frame(id=1:length(mat_data$tavg),
                         mat=mat_resid[],
                         win_tmin=win_tmin_resid[],
                         sum_tmax=sum_tmax_resid[],
                         ext_range=ann_range[])
clust_data_comp <- clust_data[complete.cases(clust_data),]
climate_k1 <- kmeans(x=clust_data_comp[,-1],centers=5,iter.max=10)
classes <- data.frame(id=clust_data_comp$id,
                      clust=fitted(climate_k1,method="classes"))
clust_data_merge <- merge(clust_data,classes,all.x=TRUE)
kmeans_clust <- ann_range
kmeans_clust[] <- clust_data_merge$clust
plot(kmeans_clust)

##Divides mean annual temperature into equal-interval bins.
class_min <- -1
mat_min <- cellStats(ann_tavg[[1]],stat="min")
mat_max <- cellStats(ann_tavg[[1]],stat="max")
mat_step= (mat_max - class_min) / 5
mat_seq <- seq(class_min,mat_max,by=mat_step)
rcl_mat <- matrix(c(mat_min-0.1,mat_seq[2],10,
                mat_seq[2],mat_seq[3],20,
                mat_seq[3],mat_seq[4],30,
                mat_seq[4],mat_seq[5],40,
                mat_seq[5],mat_seq[6]+0.1,50),
              ncol=3,byrow=TRUE)
mat_classed <- reclassify(ann_tavg[[1]],rcl=rcl_mat)
clim_classes <- mat_classed+kmeans_clust
writeRaster(clim_classes,filename="~/GIS/MORAclim_classes.tif",overwrite=TRUE)
writeRaster(mat_classed,filename="~/GIS/MORAclim_mat_classes.tif",overwrite=TRUE)
writeRaster(kmeans_clust,filename="~/GIS/MORAclim_resid_classes.tif",overwrite=TRUE)

##Estimates effect of canopy on extreme maximum temperatures.
ext_dat <- data.frame(ext_tmax=ext_tmax_resid[],
                      ext_tmin=ext_tmin_resid[],
                      tmax=sum_tmax_resid[],
                      tmin=win_tmin_resid[],
                      tavg=mat_resid[],
                      canvol=preds_90m$mcpred_canvol_81m[],
                      cair=preds_90m$mcpred_reg_coldair_3m[],
                      utmx=preds_90m$mcpred_utmx_3m[],
                      utmy=preds_90m$mcpred_utmy_3m[],
                      srad_sum=preds_90m$mcpred_srad_JJA_3m[])
ext_dat$canvol[ext_dat$canvol<0] <- 0
plot(ext_tmax~canvol,data=ext_dat,pch=".")
plot(tmax~canvol,data=ext_dat,pch=".")
plot(tmin~canvol,data=ext_dat,pch=".")
plot(ext_tmin~canvol,data=ext_dat,pch=".")

##Generalized additive model.
library(mgcv)

can_gam_tavg <- gam(tavg~s(canvol,k=4)+utmx+utmy+srad_sum+s(cair,k=4),data=ext_dat)
par(mar=c(4,4,2,2),oma=c(2,2,0,0))
plot(can_gam_tavg)

can_gam_tmax <- gam(tmax~s(canvol,k=4)+utmx+utmy+srad_sum+s(cair,k=4),data=ext_dat)
par(mar=c(4,4,2,2),oma=c(2,2,0,0))
plot(can_gam_tmax)

can_gam_ext_tmax <- gam(ext_tmax~s(canvol,k=4)+utmx+utmy+srad_sum+s(cair,k=4),data=ext_dat)
par(mar=c(4,4,2,2),oma=c(2,2,0,0))
plot(can_gam_ext_tmax)

can_gam_tmin <- gam(tmin~s(canvol,k=4)+utmx+utmy+srad_sum+s(cair,k=4),data=ext_dat)
par(mar=c(4,4,2,2),oma=c(2,2,0,0))
plot(can_gam_tmin)

can_gam_ext_tmin <- gam(ext_tmin~s(canvol,k=4)+utmx+utmy+srad_sum+s(cair,k=4),data=ext_dat)
par(mar=c(4,4,2,2),oma=c(2,2,0,0))
plot(can_gam_ext_tmin)


can_pred <- data.frame(canvol=seq(min(ext_dat$canvol,na.rm=TRUE),max(ext_dat$canvol,na.rm=TRUE),by=100),
                       srad_sum=median(ext_dat$srad_sum,na.rm=TRUE),
                       cair=median(ext_dat$cair,na.rm=TRUE),
                       utmx=median(ext_dat$utmx,na.rm=TRUE),
                       utmy=median(ext_dat$utmy,na.rm=TRUE))
can_pred$tmax_pred <- predict(can_gam_tmax,newdata=can_pred)
can_pred$ext_tmax_pred <- predict(can_gam_ext_tmax,newdata=can_pred)
can_pred$ext_tmin_pred <- predict(can_gam_ext_tmin,newdata=can_pred)
can_pred$tmin_pred <- predict(can_gam_tmin,newdata=can_pred)
can_pred$tavg_pred <- predict(can_gam_tavg,newdata=can_pred)

cair_pred <- data.frame(canvol=median(ext_dat$canvol,na.rm=TRUE),
                       srad_sum=median(ext_dat$srad_sum,na.rm=TRUE),
                       cair=seq(min(ext_dat$cair,na.rm=TRUE),max(ext_dat$cair,na.rm=TRUE),by=0.01),
                       utmx=median(ext_dat$utmx,na.rm=TRUE),
                       utmy=median(ext_dat$utmy,na.rm=TRUE))
cair_pred$tmax_pred <- predict(can_gam_tmax,newdata=cair_pred)
cair_pred$ext_tmax_pred <- predict(can_gam_ext_tmax,newdata=cair_pred)
cair_pred$ext_tmin_pred <- predict(can_gam_ext_tmin,newdata=cair_pred)
cair_pred$tmin_pred <- predict(can_gam_tmin,newdata=cair_pred)
cair_pred$tavg_pred <- predict(can_gam_tavg,newdata=cair_pred)


##Plots predictions.
pdf("Can_vol_microclim.pdf",width=5,height=6)
plot(can_pred$canvol/10000,can_pred$tmax_pred,type="l",
     xlim=c(0,80),ylim=c(-1,1.6),col=2,
     xlab="Canopy Volume (m^3 / 10000)",ylab="Temperature Anomaly (C)")
points(can_pred$canvol/10000,can_pred$ext_tmax_pred,type="l",lty=3,col=2,lwd=2)
points(can_pred$canvol/10000,can_pred$tmin_pred,type="l",lty=1,col=4,lwd=1)
points(can_pred$canvol/10000,can_pred$ext_tmin_pred,type="l",lty=3,col=4,lwd=2)
points(can_pred$canvol/10000,can_pred$tavg_pred,type="l",lwd=3,col=3)
legend("top",legend=c("Average Summer Tmax",
                           "95th Percentile Summer Tmax",
                           "Average Winter Tmin",
                           "5th Percentile Winter Tmin",
                           "Mean Annual Tavg"),
       lty=c(1,3,1,3,1),col=c(2,2,4,4,3),lwd=c(1,2,1,2,2),bty="n")
dev.off()

pdf("Cold_air_microclim.pdf",width=5,height=6)
plot(cair_pred$cair,cair_pred$tmax_pred,type="l",
     xlim=c(0,1),ylim=c(-1,1.6),col=2,
     xlab="Cold-air Index (unitless)",ylab="Temperature Anomaly (C)")
points(cair_pred$cair,cair_pred$ext_tmax_pred,type="l",lty=3,col=2,lwd=2)
points(cair_pred$cair,cair_pred$tmin_pred,type="l",lty=1,col=4,lwd=1)
points(cair_pred$cair,cair_pred$ext_tmin_pred,type="l",lty=3,col=4,lwd=2)
points(cair_pred$cair,cair_pred$tavg_pred,type="l",lty=1,col=3,lwd=3)
legend("top",legend=c("Average Summer Tmax",
                       "95th Percentile Summer Tmax",
                       "Average Winter Tmin",
                       "5th Percentile Winter Tmin",
                       "Mean Annual Tavg"),
       lty=c(1,3,1,3,1),col=c(2,2,4,4,3),lwd=c(1,2,1,2,3),bty="n")
dev.off()



