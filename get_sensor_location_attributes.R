####Script to merge sensor locations and environmental attributes with snow duration data.
####author: Ian Breckheimer
####date: 2 February 2014

####Loads required packages
library(raster)
library(maptools)
library(maps)
library(ggplot2)

##Sets the working directory.
setwd("~/code/MORA_microclimate")

##Reads in the raster attribute data
env3 <- brick("/Volumes/ib_working/GIS/env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)")
names(env3) <- env_names

##Imports aspect data.
env3$asp <- raster('/Volumes/ib_working/GIS/MORA_aspect_3m.tif')
env3$tri <- raster('/Volumes/ib_working/GIS/MORA_TRI_3m.tif')
env3$cold <- raster('/Volumes/ib_working/GIS/MORA_coldair_3m.tif')

##Reads in the coordinates of all of the microclimate sensors.
setwd("~/Dropbox/lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
sensors <- read.csv("sensor_locations_updated_10-5-2014.csv")
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

##Reads in the snow duration data and merges it with the coordinates.
attrib <- sensors_attrib@data
snow <- read.csv("Snow_cover_meta_final_2_6_2013.csv",sep=",",header=T)
snow_coords <- merge(snow,attrib,all.x=TRUE,by.x="locname_standard",by.y="combined_name")

##Computes heat load from slope and aspect data, equation from McCune and Dylan (2002)
heat_ld <- function(lat,asp,slope){
  ##Convert to radians
  lat <- lat * (2*pi) / 360
  asp <- asp * (2*pi) / 360
  slope = slope * (2*pi) / 360;
  
  ##Compute folded aspect
  a_fh <- abs(pi-abs(asp-(5*pi/4)))
  
  ##Fit parameters
  cs <- c(0.339,0.808,0.0,0.196,0,0.482)
  
  hli <- cs[1] + (cs[2]*cos(lat)*cos(slope)) + (cs[3]*cos(a_fh)*sin(slope)*sin(lat)) + (cs[4]*sin(lat)*sin(slope)) + (cs[5]*sin(a_fh)*sin(slope)) + (cs[6]*cos(a_fh)*sin(slope))
  return(hli)
}

snow_coords$heat_ld <- heat_ld(snow_coords$Lat,snow_coords$asp,snow_coords$slope)
diss_days <- strptime(as.character(snow_coords$snow_disappearance_date),format="%m/%d/%y")
snow_coords$snow_diss_DOY <- as.numeric(strftime(diss_days,format="%j"))

##Simplifies the output.
snow_output <- with(snow_coords,data.frame(study=Study,
                                           water_year=year,
                                           site_name=site_name,
                                           sensor_name=locname_standard,
                                           longitude=Long,
                                           latitude=Lat,
                                           X_UTM=X_UTM,
                                           Y_UTM=Y_UTM,
                                           elevation=elev,
                                           snow_app_date=snow_appearance_date,
                                           snow_dis_date=snow_disappearance_date,
                                           snow_cover_days=snow_cover_duration,
                                           snow_diss_DOY=snow_diss_DOY,
                                           canopy_pct=can_pct,
                                           canopy_pct_30m=canopy_pct_30m,
                                           relev_30m=relev30,
                                           relev_90m=relev90,
                                           relev_270m=relev270,
                                           tri=tri,
                                           stream_dist=stream_dist,
                                           slope=slope,
                                           aspect=asp,
                                           srad=srad,
                                           hld=heat_ld))
snow_output <- snow_output[complete.cases(snow_output),]
snow_output <- unique(snow_output)

##Quick plots to check data.
qplot(snow_cover_days,elevation,data=snow_output,facets=water_year~.,color=canopy_pct*100,xlim=c(0,365))+
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y~poly(x,2))+
  labs(x="Snow Duration (Days)",
       y="Elevation (m)",
       color="Canopy\nCover (%)",
       title="Water Year")+
  theme_bw()

##Quick plots to check data.
qplot(snow_diss_DOY,elevation,data=snow_output,facets=water_year~.,color=canopy_pct*100)+
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y~poly(x,2))+
  labs(x="Snow Dissapearance Date (Julian Day)",
       y="Elevation (m)",
       color="Canopy\nCover (%)",
       title="Water Year")+
  theme_bw()

##Exports the data to csv.
write.csv(snow_output,"microclimate_snow_coords_final_2_12_14.csv")
