##Script to extract climate variables for MORA forest stands.

library(raster)
library(maptools)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/")

##Climate Grids.
MAT <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_MAT.tif")
MAT_se <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_disp_TAVG_se.tif")
WTMIN <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_WTMIN.tif")
WTMIN_se <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_disp_WTMIN_se.tif")
DTMAX <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_STMAX.tif")
DTMAX_se <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_disp_STMAX_se.tif")
MAT_sens <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_sens_MAT.tif")
MAT_sens_se <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_sens_TAVG_se.tif")
WTMIN_sens <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_sens_WTMIN.tif")
WTMIN_sens_se <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_sens_WTMIN_se.tif")
DTMAX_sens <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_sens_STMAX.tif")
DTMAX_sens_se <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/Micro_sens_STMAX_se.tif")
PRISM_MAT <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/PRISM_MAT.tif")
PRISM_WTMIN <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/PRISM_WTMIN.tif")
PRISM_DTMAX <- raster("/Volumes/ib_working/GIS/microclim_2011_2015/PRISM_WTMIN.tif")

##Stand coordinates
stands <- read.csv("./raw/MORA STANDS_LatLongs_UTM.csv")
coordinates(stands) <- ~UTME+UTMN
crs(stands) <- crs(MAT)


clim_stack <- stack(MAT,MAT_se,
                    WTMIN,WTMIN_se,
                    DTMAX,DTMAX_se,
                    MAT_sens,MAT_sens_se,
                    WTMIN_sens,WTMIN_sens_se,
                    DTMAX_sens,DTMAX_sens_se,
                    PRISM_MAT,PRISM_WTMIN,PRISM_DTMAX)
names(clim_stack) <- c("Micro_TAVG","Micro_TAVG_se",
                       "Micro_WTMIN","Micro_WTMIN_se",
                       "Micro_DTMAX","Micro_DTMAX_se",
                       "Micro_TAVG_sens","Micro_TAVG_sens_se",
                       "Micro_WTMIN_sens","Micro_WTMIN_sens_se",
                       "Micro_DTMAX_sens","Micro_DTMAX_sens_se",
                       "PRISM_TAVG","PRISM_WTMIN","PRISM_DTMAX")

GIS_elev <- raster("/Volumes/ib_working/GIS/MORA_elev_3m.tif")
names(GIS_elev) <- "GIS_elev"

##Samples rasters at the stand coordinates.
stands_clim <- extract(clim_stack,stands,method="bilinear",sp=TRUE)
stands_clim$GIS_elev <- extract(GIS_elev,stands,method="simple",sp=FALSE,fun=mean,buffer=10)

write.csv(stands_clim,"./cleaned/stands_clim.csv",row.names=FALSE)
