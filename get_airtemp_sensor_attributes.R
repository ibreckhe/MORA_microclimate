##Script to get attributes of microclimate sensor locations.
library(raster)

sensors <- read.csv("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/sensor_locations_updated_5-11-2016.csv")

##Reads in the raster attribute data
env3 <- brick("/Volumes/ib_working/GIS/env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)")
names(env3) <- env_names

##Imports aspect and other supplementary data.
env3$asp <- raster('/Volumes/ib_working/GIS/MORA_aspect_3m.tif')
env3$tri <- raster('/Volumes/ib_working/GIS/MORA_TRI_3m.tif')
env3$cold <- raster('/Volumes/ib_working/GIS/MORA_coldair_3m.tif')
env3$srad_noc <- raster('/Volumes/ib_working/GIS/MORA_srad_yearsum_9m.tif')

elev_9m <- raster("/Volumes/ib_working/GIS/elev_NED_9m.tif")
cold_9m <- raster("/Volumes/ib_working/GIS/coldair_index_regional_9m.tif")
reg_cancov_27m <- raster("/Volumes/ib_working/GIS/reg_canopy_focal27m.tif")
reg_cancov_81m <- raster("/Volumes/ib_working/GIS/reg_canopy_focal81m.tif")
reg_cancov_243m <- raster("/Volumes/ib_working/GIS/reg_canopy_focal243m.tif")
canvol_81m <- raster('/Volumes/ib_working/GIS/mcpred_canvol_81m.tif')
cancov_81m <- raster('/Volumes/ib_working/GIS/mcpred_canpct_81m.tif')
dry_81m <- raster('/Volumes/ib_working/GIS/MORA_dry_index_81m.tif')


##Gets seasonal solar radiation.
srad_reg_DJF <- raster("/Volumes/ib_working/GIS/srad_seas_DJF.tiff")
srad_reg_MAM <- raster("/Volumes/ib_working/GIS/srad_seas_MAM.tiff")
srad_reg_JJA <- raster("/Volumes/ib_working/GIS/srad_seas_JJA.tiff")
srad_reg_SON <- raster("/Volumes/ib_working/GIS/srad_seas_SON.tiff")

##Gets UTM coordinates for the sensors.
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

sensors_attrib@data$elev <- extract(elev_9m,sensors_UTM,method="simple")
sensors_attrib@data$cair_reg_9m <- extract(cold_9m,sensors_UTM,method="simple")
sensors_attrib@data$reg_cancov_27m <- extract(reg_cancov_27m,sensors_UTM,method="simple")
sensors_attrib@data$reg_cancov_81m <- extract(reg_cancov_81m,sensors_UTM,method="simple")
sensors_attrib@data$reg_cancov_243m <- extract(reg_cancov_243m,sensors_UTM,method="simple")
sensors_attrib@data$canvol_81m <- extract(canvol_81m,sensors_UTM,method="simple")
sensors_attrib@data$cancov_81m <- extract(cancov_81m,sensors_UTM,method="simple")
sensors_attrib@data$dry_index_81m <- extract(dry_81m,sensors_UTM,method="simple")
sensors_attrib@data$srad_seas_DJF <- extract(srad_reg_DJF,sensors_UTM,method="simple")
sensors_attrib@data$srad_seas_MAM <- extract(srad_reg_MAM,sensors_UTM,method="simple")
sensors_attrib@data$srad_seas_JJA <- extract(srad_reg_JJA,sensors_UTM,method="simple")
sensors_attrib@data$srad_seas_SON <- extract(srad_reg_SON,sensors_UTM,method="simple")

##Computes the average, solar radiation, canopy height and percentage within 30m
sensors_attrib@data$can_ht_30m <- extract(env3$can_ht,sensors_UTM,buffer=30,fun=mean)
sensors_attrib@data$canopy_pct_30m <- extract(env3$can_pct,sensors_UTM,buffer=30,fun=mean)
sensors_attrib@data$srad_can_30m <- extract(env3$srad,sensors_UTM,buffer=30,fun=mean)

##Subsets air tempertature sites.
sensors_airtemp <- subset(sensors_attrib,sub_site2 %in% c("A1","A2","A-1","A-2"))
sensors_airtemp$type <- "DataLogger"
sensors_airtemp$type[sensors_airtemp$study %in% c("MesoWest","NPSClimate")] <- "WeatherStation"
write.csv(sensors_airtemp,"/Users/ian/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/airtemp_sensor_metadata_2016.csv",row.names=FALSE)
