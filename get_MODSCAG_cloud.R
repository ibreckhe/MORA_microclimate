##Script to extract MODSCAG cloud fraction for each day in the microclimate dataset.
##Author: Ian Breckheimer
##October 30th, 2014

##Sets up workspace.
library(R.matlab)
library(raster)
library(rgdal)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Nicoleta")

##Reads in data (takes a while).
mat_dat <- readMat("MSG_data_cube_2000_2013_raw.mat")
names(mat_dat)

##Reads into a raster brick object.
modscag_brick <- brick(mat_dat$fSCA,xmn=min(mat_dat$scene.X)-231.65,xmx=max(mat_dat$scene.X)+231.65,
                       ymn=min(mat_dat$scene.Y)-231.5,ymx=max(mat_dat$scene.Y)+231.5)
projection(modscag_brick) <- "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
writeRaster(modscag_brick,filename="MSG_data_cube_2000_2013.grd",datatype="INT1U",overwrite=TRUE)

##Reads in data and gets rid of values other than clouds and snow.
msg <- brick("MSG_data_cube_2000_2013.grd")
days <- format(seq(from=as.Date("2000-10-01"),to=as.Date("2014-10-01"),by="day"),format="%b.%d.%Y")
names(msg) <- days
rcl_df <- data.frame(id=c(0:100,230,235,250,253,255),v=c(rep(0,101),NA,100,100,NA,NA))
mst_rcl <- subs(msg[[209:309]],rcl_df,filename="MSG_cloud.grd",datatype="INT2U",overwrite=TRUE)
mst_stats <- cellStats(mst_rcl,stat="mean")
plot(1:length(mst_stats),mst_stats,type="l")
