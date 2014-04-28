##Script to predict snow disapearance dates of meadowatch sites.
##Author:Ian Breckheimer
##Date: 14 February 2014

##Loads required packages
library(raster)
library(ggmap)

##Loads data
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/MeadoWatch")
msites <- read.csv("Locations of MeadoWatch Sites.csv")
msites
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/flickr_inat")
fsites <- read.csv("Flickr_iNat_all_through_2013.csv")
head(fsites)
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/Flickr/api_processing")
asites <- read.csv("MORA_flickr_all_2001_2013.csv")
head(asites)
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/Flickr/")
f2013_sites <- read.csv("MORA_flickr_classified_all_2013.csv")
head(f2013_sites)

##Puts them on a map to check coordinates
ymx <- max(msites$latitude_dd)
ymn <- min(msites$latitude_dd)
xmx <- max(msites$longitude_dd)
xmn <- min(msites$longitude_dd)
bbox <- c(xmn,ymn,xmx,ymx)
map <- get_map(bbox,zoom=15,maptype="terrain",source="google",color="bw")
ggmap(map)+geom_point(aes(x = msites$longitude_dd,y = msites$latitude_dd,color=msites$station))

##Loads SDD predictions
setwd("~/GIS")
snowmaps <- c("MORA_snow_pred_2009_3m.tif",
              "MORA_snow_pred_2010_3m.tif",
              "MORA_snow_pred_2011_3m.tif",
              "MORA_snow_pred_2012_3m.tif",
              "MORA_snow_pred_2013_3m.tif")
snowstack <- stack(snowmaps)
names(snowstack) <- paste(rep(c("estimate","lwr","upr"),5),sort(rep(2009:2013,3)),sep="_")

##Brings in the 3m environmental data
setwd("/Users/ian/GIS/")
env3 <- brick("env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)")
names(env3) <- env_names

##Converts locations to a spatial points data frame.
coordinates(msites) <- ~UTM_X+UTM_Y
proj4string(msites) <- "+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
coordinates(fsites) <- ~UTME+UTMN
proj4string(fsites) <- "+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
fsites@data$UTM_E <- coordinates(fsites)[,1]
fsites@data$UTM_N <- coordinates(fsites)[,2]
coordinates(asites) <- ~longitude+latitude
proj4string(asites) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"
asites <- spTransform(asites,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
asites@data$UTM_E <- coordinates(asites)[,1]
asites@data$UTM_N <- coordinates(asites)[,2]

##Extracts melt dates with standard errors for all of the sites.
fsites_snow <- extract(snowstack,fsites,method="simple",sp=TRUE)
msites_snow <- extract(snowstack,msites,method="simple",sp=TRUE)
asites_snow <- extract(snowstack,asites,method="simple",sp=TRUE)


##Samples the other covariate data at the points.
fsites_covar <- extract(env3,fsites_snow,method="simple",sp=TRUE)
msites_covar <- extract(env3,msites_snow,method="simple",sp=TRUE)
asites_covar <- extract(env3,asites_snow,method="simple",sp=TRUE)

##Extracts canopy cover at a 30m radius
#fsites_covar@data$can_pct_30m <- extract(env3[["can_pct"]],fsites_snow,buffer=30,fun=mean)
#msites_covar@data$can_pct_30m <- extract(env3[["can_pct"]],msites_snow,buffer=30,fun=mean)
#asites_covar@data$can_pct_30m <- extract(env3[["can_pct"]],asites_snow,buffer=30,fun=mean)

##Constructs a days since snow field
fsites_2009 <- fsites_covar[fsites_covar$Year==2009,]
fsites_2010 <- fsites_covar[fsites_covar$Year==2010,]
fsites_2011 <- fsites_covar[fsites_covar$Year==2011,]
fsites_2012 <- fsites_covar[fsites_covar$Year==2012,]
fsites_2013 <- fsites_covar[fsites_covar$Year==2013,]

fsites_2009$dss <- fsites_2009$Day.of.Year - fsites_2009$estimate_2009
fsites_2009$dss_upr <- fsites_2009$Day.of.Year - fsites_2009$lwr_2009
fsites_2009$dss_lwr <- fsites_2009$Day.of.Year - fsites_2009$upr_2009

fsites_2010$dss <- fsites_2010$Day.of.Year - fsites_2010$estimate_2010
fsites_2010$dss_upr <- fsites_2010$Day.of.Year - fsites_2010$lwr_2010
fsites_2010$dss_lwr <- fsites_2010$Day.of.Year - fsites_2010$upr_2010

fsites_2011$dss <- fsites_2011$Day.of.Year - fsites_2011$estimate_2011
fsites_2011$dss_upr<- fsites_2011$Day.of.Year - fsites_2011$lwr_2011
fsites_2011$dss_lwr <- fsites_2011$Day.of.Year - fsites_2011$upr_2011

fsites_2012$dss <- fsites_2012$Day.of.Year - fsites_2012$estimate_2012
fsites_2012$dss_upr <- fsites_2012$Day.of.Year - fsites_2012$lwr_2012
fsites_2012$dss_lwr <- fsites_2012$Day.of.Year - fsites_2012$upr_2012

fsites_2013$dss <- fsites_2013$Day.of.Year - fsites_2013$estimate_2013
fsites_2013$dss_upr <- fsites_2013$Day.of.Year - fsites_2013$lwr_2013
fsites_2013$dss_lwr <- fsites_2013$Day.of.Year - fsites_2013$upr_2013

fsites_dss <- rbind(fsites_2009,fsites_2010,fsites_2011,fsites_2012,fsites_2013)
fsites_dss@data[,44] <- as.numeric(fsites_dss@data[,44])


##Constructs a days since snow field for all of the sites.
asites_covar$datetaken <- as.Date(as.character(asites_covar$datetaken))
asites_covar$Year <- as.numeric(strftime(asites_covar$datetaken,format="%Y"))
asites_covar$Day.of.Year <- as.numeric(strftime(asites_covar$datetaken,format="%j"))

asites_2009 <- asites_covar[asites_covar$Year==2009,]
asites_2010 <- asites_covar[asites_covar$Year==2010,]
asites_2011 <- asites_covar[asites_covar$Year==2011,]
asites_2012 <- asites_covar[asites_covar$Year==2012,]
asites_2013 <- asites_covar[asites_covar$Year==2013,]

asites_2009$dss <- asites_2009$Day.of.Year - asites_2009$estimate_2009
asites_2009$dss_upr <- asites_2009$Day.of.Year - asites_2009$lwr_2009
asites_2009$dss_lwr <- asites_2009$Day.of.Year - asites_2009$upr_2009

asites_2010$dss <- asites_2010$Day.of.Year - asites_2010$estimate_2010
asites_2010$dss_upr <- asites_2010$Day.of.Year - asites_2010$lwr_2010
asites_2010$dss_lwr <- asites_2010$Day.of.Year - asites_2010$upr_2010

asites_2011$dss <- asites_2011$Day.of.Year - asites_2011$estimate_2011
asites_2011$dss_upr<- asites_2011$Day.of.Year - asites_2011$lwr_2011
asites_2011$dss_lwr <- asites_2011$Day.of.Year - asites_2011$upr_2011

asites_2012$dss <- asites_2012$Day.of.Year - asites_2012$estimate_2012
asites_2012$dss_upr <- asites_2012$Day.of.Year - asites_2012$lwr_2012
asites_2012$dss_lwr <- asites_2012$Day.of.Year - asites_2012$upr_2012

asites_2013$dss <- asites_2013$Day.of.Year - asites_2013$estimate_2013
asites_2013$dss_upr <- asites_2013$Day.of.Year - asites_2013$lwr_2013
asites_2013$dss_lwr <- asites_2013$Day.of.Year - asites_2013$upr_2013

asites_dss <- rbind(asites_2009,asites_2010,asites_2011,asites_2012,asites_2013)

##Filters flickr and iNaturalist data based on covariates.
fsites_clean <- subset(fsites_dss,elev > 800 & 
                                    elev < 2500 & 
                                    snow_melt != 2)
##Filters background data based on covariates.
asites_clean <- subset(asites_dss,elev > 400)

##Merges the flickr databases together by common columns.
asites_merge <- data.frame(datetaken=asites_clean@data$datetaken,
                           DOY=asites_clean@data$Day.of.Year,
                           Year=asites_clean@data$Year,
                           Species="NF",
                           Phenological.Stage="NF",
                           Flowering=0,
                           Title=asites_clean@data$title,
                           ID=asites_clean@data$id,
                           Owner=asites_clean@data$owner,
                           ID.Person=NA,
                           Source="Flickr")
asites_merge <- cbind(asites_merge,asites_clean@data[,c(30:60,63:65)])
fsites_merge <- fsites_clean@data
fsites_merge$Date <- strptime(fsites_clean@data$Date,format="%m/%d/%y")
colnames(fsites_merge) <- colnames(asites_merge)

asites_all <- rbind(asites_merge,fsites_merge)

##Writes to csv file
setwd("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/cleaned")
write.csv(msites_covar@data,"meadowatch_locations_snow_covar.csv",row.names=FALSE)
write.csv(fsites_clean@data,"flickr_inat_locations_snow_filtered_dss.csv",row.names=FALSE)
write.csv(asites_clean@data,"flickr_locations_all_2001_2013_snow_filtered.csv",row.names=FALSE)
write.csv(asites_all,"flickr_inat_all_classed_2001_2013.csv",row.names=FALSE)

