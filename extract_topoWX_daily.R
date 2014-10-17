## Script to extract daily time-series from topoWX 800m data.
## Author: Ian Breckheimer.

## Sets up workspace.
library(raster)
library(ncdf4)

## Creates a bounding box around Mt. Rainier.
ran_ext <- extent(matrix(c(-122.172553,-120.914619,46.487004,47.212064),
                         ncol=2,byrow=TRUE))

## Opens the tmax ncdf file.
tmax_nc <- nc_open("~/Downloads/topoWX_tmax_2004_2012.nc")
tmax_var <- "tmax"

## Gets coordinates of cell corners
lon <- ncvar_get(tmax_nc,"lon")
xmin <- min(lon) - (0.5*(lon[2]-lon[1]))
xmax <- max(lon) + (0.5*(lon[2]-lon[1]))
lat <- ncvar_get(tmax_nc,"lat")
ymin <- min(lat) - (0.5*(lat[2]-lat[1]))
ymax <- max(lat) + (0.5*(lat[2]-lat[1]))
tmp_array <- ncvar_get(tmax_nc,"tmax")
tmp_raster <- brick(tmp_array,transpose=TRUE,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,
                    crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
crop_raster <- crop(tmp_raster,y=ran_ext,filename="~/Desktop/topowx_tmax_2004-2012.grd")

## Extracts time-series of daily averages.
tmax_daily <- cellStats(crop_raster,'mean')

## Opens the tmin ncdf file.
tmin_nc <- nc_open("~/Downloads/topoWX_tmin_2004_2012.nc")
tmin_var <- "tmin"

## Gets coordinates of cell corners
lon <- ncvar_get(tmin_nc,"lon")
xmin <- min(lon) - (0.5*(lon[2]-lon[1]))
xmax <- max(lon) + (0.5*(lon[2]-lon[1]))
lat <- ncvar_get(tmax_nc,"lat")
ymin <- min(lat) - (0.5*(lat[2]-lat[1]))
ymax <- max(lat) + (0.5*(lat[2]-lat[1]))
tmp_array <- ncvar_get(tmin_nc,"tmin")
tmp_raster <- brick(tmp_array,transpose=TRUE,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,
                    crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
crop_raster <- crop(tmp_raster,y=ran_ext,filename="~/Desktop/topowx_tmin_2004-2012.grd")

## Extracts time-series of daily averages.
tmin_daily <- cellStats(crop_raster,'mean')

## Calculates a daily tavg series, and writes all the data to disk.
tavg_daily <- (tmax_daily + tmin_daily) / 2
days <- seq(from=as.Date("2004-01-01"),to=as.Date("2012-12-31"),by="day")
topowx_daily <- data.frame(DATE=days,TMAX=tmax_daily,TMIN=tmin_daily,TAVG=tavg_daily)
write.csv(topowx_daily,"topowx_daily_2004_2012.csv",row.names=FALSE)


