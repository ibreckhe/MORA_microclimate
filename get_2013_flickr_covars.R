##Script to extract predicted snow melt and other spatial covariates for flickr points.
##Ian Breckheimer
##4-9-2014

##Sets up the workspace
library(raster)
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/Flickr/")

##Loads flickr data and converts to spatial coordinates.
f2013 <- read.csv("MORA_flickr_classified_all_2013.csv")
coordinates(f2013) <- ~UTME+UTMN
crs(f2013) <- "+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
f2013@data$UTM_X <- coordinates(f2013)[,1]
f2013@data$UTM_Y <- coordinates(f2013)[,2]
head(f2013)

##Loads GIS data
env3 <- brick("~/GIS/env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
##Imports aspect data and computes relative elevation.
env3$asp <- raster('~/GIS/MORA_aspect_3m.tif')
##Computes relative elevation (elevation of cell - average within 30m)
gf <- focalWeight(env3[[4]], 30, "circle")
env3$relev_30m <- env3[[4]] - focal(env3[[4]],gf)
##adds predicted snow melt
env3$predsnow2013 <- raster("~/GIS/MORA_snow_pred_2013_lme.tif")

##Cleans up names
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist","asp","relev_30m","snow_diss_2013")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Summer Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)","Aspect (deg)",
                "Relative Elevation (30m)","Predicted SDD 2013")
names(env3) <- env_names

##Samples all of the covariates at the sensor locations.
f2013_attrib <- extract(env3,f2013,method='simple',sp=T)

##Converts coordinates to lat-long.
f2013_attrib_latlong <- spTransform(f2013_attrib,CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
f2013_attrib_latlong@data$long <- coordinates(f2013_attrib_latlong)[,1]
f2013_attrib_latlong@data$lat <- coordinates(f2013_attrib_latlong)[,2]

##Saves the flickr data to disk.
write.csv(f2013_attrib_latlong@data,"MORA_flickr_classified_all_2013_covars.csv",row.names=FALSE)

