####Script to append GIS attributes such as slope and elevation to the file with site geolocations.
####author: Ian Breckheimer
####date: 14 January 2014

####Loads required packages
library(raster)
library(maptools)
library(maps)

##Sets the working directory.
setwd("~/code/MORA_microclimate")

##Reads in the raster attribute data
env3 <- brick("~/GIS/env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)")
names(env3) <- env_names

##Reads in the coordinates.
sensors <- read.csv("~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/sensor_locations.csv")
coordinates(sensors) <- ~Long+Lat
sensors@data$Long <- coordinates(sensors)[,1]
sensors@data$Lat <- coordinates(sensors)[,2]
crs(sensors) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"

##Transforms the coordinates to UTM and appends coordinates to the data frame.
sensors_UTM <- spTransform(sensors,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
sensors_UTM@data$X_UTM <- coordinates(sensors_UTM)[,1]
sensors_UTM@data$Y_UTM <- coordinates(sensors_UTM)[,2]

##Samples the raster brick at the sensor locations.
sensors_attrib <- extract(env3,sensors_UTM,method='simple',sp=T)

##Computes relative elevation (elevation of cell - average within 90m)
relev30 <- extract(env3$elev,sensors_UTM,buffer=30,fun=mean)
relev90 <- extract(env3$elev,sensors_UTM,buffer=90,fun=mean)
relev270 <- extract(env3$elev,sensors_UTM,buffer=270,fun=mean)
sensors_attrib@data$relev30 <- sensors_attrib@data$elev - relev30
sensors_attrib@data$relev90 <- sensors_attrib@data$elev - relev90
sensors_attrib@data$relev270 <- sensors_attrib@data$elev - relev270

##Computes the average canopy height and percentage within 30m
sensors_attrib@data$can_ht_30m <- extract(env3$can_ht,sensors_UTM,buffer=30,fun=mean)
sensors_attrib@data$canopy_pct_30m <- extract(env3$can_pct,sensors_UTM,buffer=30,fun=mean)

##Exports the data to csv.
write.csv(sensors_attrib@data,
          "~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/sensor_locations_attrib.csv")
