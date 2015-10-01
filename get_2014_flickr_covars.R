##Script to extract predicted snow melt and other spatial covariates for flickr points.
##Ian Breckheimer
##4-4-2014

##Sets up the workspace
library(raster)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/Flickr/")

##Loads flickr data and converts to spatial coordinates.
f <- read.csv("MORA_flickr_metadata_all_2009_2014.csv")
coordinates(f) <- ~UTME+UTMN
crs(f) <- "+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
f@data$UTM_X <- coordinates(f)[,1]
f@data$UTM_Y <- coordinates(f)[,2]
head(f)

##Loads GIS data
env3 <- brick("/Volumes/ib_working/GIS/env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
##Imports aspect data and computes relative elevation.
env3$asp <- raster('/Volumes/ib_working/GIS/MORA_aspect_3m.tif')
##Computes relative elevation (elevation of cell - average within 30m)
gf <- focalWeight(env3[[4]], 30, "circle")
env3$relev_30m <- env3[[4]] - focal(env3[[4]],gf)

##adds predicted snow melt
#env3$predsnow2013 <- raster("~/GIS/MORA_snow_pred_2013_lme.tif")

##Cleans up names
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist","asp","relev_30m")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Summer Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)","Aspect (deg)",
                "Relative Elevation (30m)")
names(env3) <- env_names

##Samples all of the covariates at the sensor locations.
f_attrib <- extract(env3,f,method='simple',sp=T)

##Discards points without any metadata.
f_attrib_park <- subset(f_attrib, !is.na(elev))

##Saves the flickr data to disk.
write.csv(f_attrib_park@data,"MORA_flickr_metdata_all_2009_2014_covars.csv",row.names=FALSE)

