####Script to compare MORAclim predictions to PRISM and TopoWX data.

library(ncdf4)
library(raster)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/")

##Loads topowx netcdf file.
tmax_brick <- brick("./raw/MORA_topoWX_2009_2015.nc",varname="tmax")
tmin_brick <- brick("./raw/MORA_topoWX_2009_2015.nc",varname="tmin")
elev <- raster("/volumes/ib_working/GIS/prism/ra_clim/data/NED_dem_30sec.img")
elev_clip <- crop(elev,tmax_brick)
dayseq <- seq(as.Date("2008-12-31"),as.Date("2015-12-30"),by="day")
mc_dayseq <- seq(as.Date("2009-09-01"),as.Date("2015-08-30"),by="day")
mc_months <- as.numeric(format(MORAclim_dayseq,format="%m"))
mc_djf_dayseq <- mc_dayseq[which(mc_months %in% c(12,1,2))]
mc_mam_dayseq <- mc_dayseq[which(mc_months %in% c(3,4,5))]
mc_jja_dayseq <- mc_dayseq[which(mc_months %in% c(6,7,8))]
mc_son_dayseq <- mc_dayseq[which(mc_months %in% c(9,10,11))]


mc_indices <- which(dayseq %in% mc_dayseq)
mc_djf_indices <- which(dayseq %in% mc_djf_dayseq)
mc_mam_indices <- which(dayseq %in% mc_mam_dayseq)
mc_jja_indices <- which(dayseq %in% mc_jja_dayseq)
mc_son_indices <- which(dayseq %in% mc_son_dayseq)

##Drops days outside the time window of interest
tmax_2009_2015 <- tmax_brick[[mc_indices]]
tmin_2009_2015 <- tmin_brick[[mc_indices]]

##Computes mean temperature for the time period.
tavg_2009_2015 <- (tmax_2009_2015 + tmin_2009_2015)/2
mat_2009_2015 <- calc(tavg_2009_2015,fun=mean)

##Drops days outside the time window of interest
tmax_djf <- tmax_brick[[mc_djf_indices]]
tmin_djf <- tmin_brick[[mc_djf_indices]]
tmax_mam <- tmax_brick[[mc_mam_indices]]
tmin_mam <- tmin_brick[[mc_mam_indices]]
tmax_jja <- tmax_brick[[mc_jja_indices]]
tmin_jja <- tmin_brick[[mc_jja_indices]]
tmax_son <- tmax_brick[[mc_son_indices]]
tmin_son <- tmin_brick[[mc_son_indices]]

##Computes mean temperature for the time period.
tavg_2009_2015 <- (tmax_2009_2015 + tmin_2009_2015)/2
tavg_djf <- (tmax_djf + tmin_djf)/2
tavg_mam <- (tmax_mam + tmin_mam)/2
tavg_jja <- (tmax_jja + tmin_jja)/2
tavg_son <- (tmax_son + tmin_son)/2


mat_2009_2015 <- calc(tavg_2009_2015,fun=mean)
tavg_djf_mean <- calc(tavg_djf,fun=mean)
tavg_mam_mean <- calc(tavg_mam,fun=mean)
tavg_jja_mean <- calc(tavg_jja,fun=mean)
tavg_son_mean <- calc(tavg_son,fun=mean)

tmax_djf_mean <- calc(tmax_djf,fun=mean)
tmax_mam_mean <- calc(tmax_mam,fun=mean)
tmax_jja_mean <- calc(tmax_jja,fun=mean)
tmax_son_mean <- calc(tmax_son,fun=mean)

tmin_djf_mean <- calc(tmin_djf,fun=mean)
tmin_mam_mean <- calc(tmin_mam,fun=mean)
tmin_jja_mean <- calc(tmin_jja,fun=mean)
tmin_son_mean <- calc(tmin_son,fun=mean)

##Puts everything in a raster brick.
tx_all <- brick(mat_2009_2015,
                tavg_djf_mean,tavg_mam_mean,tavg_jja_mean,tavg_son_mean,
                tmin_djf_mean,tmin_mam_mean,tmin_jja_mean,tmin_son_mean,
                tmax_djf_mean,tmax_mam_mean,tmax_jja_mean,tmax_son_mean)

##Brings in MORAclim estimates.
mc_ann_tavg <- brick("/Volumes/ib_working/MORAclim_summaries/ann_tavg_sum.tif")[[1]]
mc_djf_tavg <- brick("/Volumes/ib_working/MORAclim_summaries/djf_tavg_sum.tif")[[1]]
mc_mam_tavg <- brick("/Volumes/ib_working/MORAclim_summaries/mam_tavg_sum.tif")[[1]]
mc_jja_tavg <- brick("/Volumes/ib_working/MORAclim_summaries/jja_tavg_sum.tif")[[1]]
mc_son_tavg <- brick("/Volumes/ib_working/MORAclim_summaries/son_tavg_sum.tif")[[1]]

mc_djf_tmin <- brick("/Volumes/ib_working/MORAclim_summaries/djf_tmin_sum.tif")[[1]]
mc_mam_tmin <- brick("/Volumes/ib_working/MORAclim_summaries/mam_tmin_sum.tif")[[1]]
mc_jja_tmin <- brick("/Volumes/ib_working/MORAclim_summaries/jja_tmin_sum.tif")[[1]]
mc_son_tmin <- brick("/Volumes/ib_working/MORAclim_summaries/son_tmin_sum.tif")[[1]]

mc_djf_tmax <- brick("/Volumes/ib_working/MORAclim_summaries/djf_tmax_sum.tif")[[1]]
mc_mam_tmax <- brick("/Volumes/ib_working/MORAclim_summaries/mam_tmax_sum.tif")[[1]]
mc_jja_tmax <- brick("/Volumes/ib_working/MORAclim_summaries/jja_tmax_sum.tif")[[1]]
mc_son_tmax <- brick("/Volumes/ib_working/MORAclim_summaries/son_tmax_sum.tif")[[1]]

mc_all <- brick(mc_ann_tavg,
                mc_djf_tavg,mc_mam_tavg,mc_jja_tavg,mc_son_tavg,
                mc_djf_tmin,mc_mam_tmin,mc_jja_tmin,mc_son_tmin,
                mc_djf_tmax,mc_mam_tmax,mc_jja_tmax,mc_son_tmax)

tx_all_hr <- projectRaster(tx_all,mc_all,
                              method="bilinear")

##Compares the two products.
mc_tx_diff <- mc_all - tx_all_hr
mora_bound <- readShapePoly("/Volumes/ib_working/GIS/MORA_boundary.shp")

par(mar=c(0,0,0,0),oma=c(2,2,0,0))
pdf("./results/MORAclim_topowx_diff_seas.pdf",width=8,height=6)
resid_pal <- colorRampPalette(colors=c("darkblue","white","red"))
plot(mc_tx_diff[[6:13]],col=resid_pal(n=length(resid_seq)),
     main="",zlim=c(-6,6),nc=4,box=FALSE,legend=FALSE,
     axes=TRUE,interpolate=FALSE,addfun=plot(mora_bound,lwd=2,add=TRUE))
scalebar(10000,xy=c(592000,5175000),type="line",adj=c(0.5,1.5),lwd=3,label=c("10km"))
plot(mc_tx_diff[[6:13]], col=resid_pal(n=length(resid_seq)),
     zlim=c(-5,5),legend.only = TRUE)
dev.off()


