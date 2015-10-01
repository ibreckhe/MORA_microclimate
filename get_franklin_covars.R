## Loads required packages.
library(raster)
library(rgdal)

## Loads and formats data.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
fsites <- read.csv("Franklindata_UTM.csv")
coordinates(fsites) <- ~UTME+UTMN
proj4string(fsites) <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

##Loads and stacks rasters
setwd("/Volumes/ib_working/GIS/")
elev <- raster("MORA_elev_3m.tif")
elev_9m <- raster("MORA_elev_focal9m.tif")
elev_81m <- raster("MORA_elev_focal81m.tif")
elev_243m <- raster("MORA_elev__focal_243m.tif")
elev_729m <- raster("MORA_elev__focal_729m.tif")
cvol_81m <- raster("MORA_can_vol_81m.tif")
ccov_81m <- raster("MORA_can_pct_focal81m.tif")
cold_ind_9m <- raster("MORA_coldair_index.tif")
srad_sum_9m <- raster("MORA_srad_yearsum_9m.tif")
topoidx_9m <- raster("MORA_topoidx_9m_log.tif")
topoidx_81m <- raster("MORA_topoidx_81m_log.tif")
dry_idx_9m <- raster("MORA_dry_index_9m_precip.tif")
dry_idx_81m <- raster("MORA_dry_index_81m_precip.tif")
snow_2013 <- raster("MORA_snow_pred_2013_3m_parkbound_est.tiff")
rstack <- stack(elev,elev_9m,elev_81m,elev_243m,elev_729m,
                cvol_81m,ccov_81m,cold_ind_9m,srad_sum_9m,
                topoidx_9m,topoidx_81m,dry_idx_9m,dry_idx_81m,snow_2013)

##Brings in the resampled PRISM data.
setwd("/Volumes/ib_working/GIS/prism/ra_clim/output/")
temp_monthly <- brick("tavgs_multiband_highres.img")
ppt_monthly <- brick("ppts_multiband_highres.img")
mat <- raster("MAT_c.img")
ppt <- raster("PPT_mm.img")
pstack <- stack(temp_monthly,ppt_monthly,mat,ppt)



##Samples rasters at plot locations.
fsites_ext <- extract(rstack,y=fsites,method="simple",sp=TRUE)
fsites_latlong <- spTransform(fsites,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
latlong_coords <- coordinates(fsites_latlong)
fsites_ext$lat <- latlong_coords[,2]
fsites_ext$lon <- latlong_coords[,1]
fsites_ext$utmx <- coordinates(fsites)[,1]
fsites_ext$utmy <- coordinates(fsites)[,2]
fsites_prsm <- extract(pstack,fsites_ext,method="simple",sp=TRUE)

##Merges with climateWNA data.
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/")
wna <- read.csv("franklin_latlong_only_climateWNA.csv")
wna$SITE <- paste(wna$ID2,wna$ID1)
fsites_prsm@data$SITE <- paste(fsites_prsm@data$SITE_NAME,fsites_prsm@data$PLOT_NO)
fsites_all_covar <- merge(fsites_prsm@data,wna)
write.csv(fsites_all_covar,"franklindata_utm_covars.csv",row.names=FALSE)
