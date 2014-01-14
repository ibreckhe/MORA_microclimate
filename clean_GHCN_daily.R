##Script to clean climate data downloaded from the global historical climatology network.
##Author: Ian Breckheimer
##15 December 2013

library(devtools)
library(rnoaa)
library(reshape2)

##Sets the API key for downloading NOAA data.
options(noaakey="NjuEQnGAERfxLrZNIiHqMUJsSeqiPuRR")

##Sets the working directory.
setwd("~/code/mora_climate")

# ##All (current and historical) stations in the vicinity of Mt. Rainier.
# stationIDs <- c("GHCND:USC00450285",
                # "GHCND:USC00450945",
                # "GHCND:USC00450969",
                # "GHCND:USS0021C38S",
                # "GHCND:USS0021B63S",
                # "GHCND:USC00451113",
                # "GHCND:USC00451110",
                # "GHCND:USS0021B13S",
                # "GHCND:USC00452493",
                # "GHCND:USC00452722",
                # "GHCND:USC00453177",
                # "GHCND:USC00454286",
                # "GHCND:USC00453219",
                # "GHCND:USR0000WGRN",
                # "GHCND:USC00453357",
                # "GHCND:USR0000WHAG",
                # "GHCND:USS0021B62S",
                # "GHCND:USC00454620",
                # "GHCND:USC00454619",
                # "GHCND:USC00456894",
                # "GHCND:USC00454764",
                # "GHCND:US1WALW0015",
                # "GHCND:USC00455425",
                # "GHCND:USS0021C17S",
                # "GHCND:USS0021C40S",
                # "GHCND:USC00455704",
                # "GHCND:USC00456201",
                # "GHCND:USC00456265",
                # "GHCND:USC00456262",
                # "GHCND:USS0021C35S",
                # "GHCND:USC00456381",
                # "GHCND:USS0021C33S",
                # "GHCND:USC00456892",
                # "GHCND:USC00456896",
                # "GHCND:USC00456898",
                # "GHCND:USC00456900",
                # "GHCND:USC00456909",
                # "GHCND:USC00458653",
                # "GHCND:USS0021C28S",
                # "GHCND:USC00459170",
                # "GHCND:USC00456385",
                # "GHCND:USC00459171")

####Brings in and cleans the historical station data.####
statdat <- read.csv("MORA_GHCN_all.csv")

##Cleans the TMAX data.
tmax_dat <- statdat[,c("STATION","LATITUDE","LONGITUDE","DATE","TMAX")]
tmax_dat$DATE <- strptime(as.character(tmax_dat$DATE),format="%Y%m%d")

##Creates a continuous sequence of days.
start <- min(tmax_dat$DATE)
end <- max(tmax_dat$DATE)
days <- data.frame(date=seq(start,end,by="day"))
days$CHARDATE <- substr(as.character(days$date),start=1,stop=10)

##Examines the data.
scatter.smooth(as.Date(tmax_dat$DATE),tmax_dat$TMAX,pch=20,cex=0.2,
               col=rgb(0.5,0.5,0.5,0.1),main="TMAX",lpars=list(col="red"))

##Eliminates impossible values
tmax_dat$TMAX[which(tmax_dat$TMAX > 500)] <- NA
tmax_dat$TMAX[which(tmax_dat$TMAX < -300)] <- NA

##Casts the tmax data to "wide" format.
tmax_wide<- dcast(tmax_dat,DATE~STATION,value.var="TMAX")
tmax_wide$CHARDATE <- as.character(tmax_wide$DATE)
rownames(tmax_wide) <- tmax_wide$CHARDATE
dim(tmax_wide)

##Adds in any missing days.
tmax_wide_all <- merge(days,tmax_wide,by="CHARDATE",all=T)

##Computes a median across all stations for each day.
tmax_median <- apply(tmax_wide_all[,4:dim(tmax_wide_all)[2]],MARGIN=1,FUN=median,na.rm=T)

##Checks the data.
plot(tmax_median,type="l")

##Eliminates all values more than 20C from the median.
for (i in 4:(dim(tmax_wide_all)[2])){
	diff <- abs(tmax_wide_all[,i] - as.numeric(tmax_median))
	out <- which(diff>200)
	tmax_wide_all[out,i] <- NA
}

####Cleans the TMIN data.####
tmin_dat <- statdat[,c("STATION","LATITUDE","LONGITUDE","DATE","TMIN")]
tmin_dat$DATE <- strptime(as.character(tmin_dat$DATE),format="%Y%m%d")

##Examines the data.
# scatter.smooth(as.Date(tmin_dat$DATE),tmin_dat$TMIN,pch=20,cex=0.2,
               # col=rgb(0.5,0.5,0.5,0.5),main="TMIN",lpars=list(col="red"))

##Eliminates impossible values
tmin_dat$TMIN[which(tmin_dat$TMIN > 250)] <- NA
tmin_dat$TMIN[which(tmin_dat$TMIN < -300)] <- NA

##Casts the tmax data to "wide" format.
tmin_wide<- dcast(tmin_dat,DATE~STATION,value.var="TMIN")
tmin_wide$CHARDATE <- as.character(tmin_wide$DATE)
rownames(tmin_wide) <- tmin_wide$CHARDATE
dim(tmin_wide)

##Adds in any missing days.
tmin_wide_all <- merge(days,tmin_wide,by="CHARDATE",all=T)

##Computes a median across all stations for each day.
tmin_median <- apply(tmin_wide_all[,4:dim(tmin_wide_all)[2]],MARGIN=1,FUN=median,na.rm=T)

##Eliminates all values more than 20C from the median.
for (i in 4:dim(tmin_wide_all)[2]){
	diff <- abs(tmin_wide_all[,i] - tmin_median)
	out <- which(diff>200)
	tmin_wide_all[out,i] <- NA
}

####Cleans the PRCP data.####
prcp_dat <- statdat[,c("STATION","LATITUDE","LONGITUDE","DATE","PRCP")]
prcp_dat$DATE <- strptime(as.character(prcp_dat$DATE),format="%Y%m%d")

##Examines the data.
# scatter.smooth(as.Date(prcp_dat$DATE),prcp_dat$PRCP,pch=20,cex=0.2,
               # col=rgb(0.5,0.5,0.5,0.5),main="PRCP",lpars=list(col="red"))

##Eliminates impossible values
prcp_dat$PRCP[which(prcp_dat$PRCP > 2000)] <- NA
prcp_dat$PRCP[which(prcp_dat$PRCP < 0)] <- NA

##Casts the prcp data to "wide" format.
prcp_wide<- dcast(prcp_dat,DATE~STATION,value.var="PRCP")
prcp_wide$CHARDATE <- as.character(prcp_wide$DATE)
dim(prcp_wide)

##Adds in any missing days.
prcp_wide_all <- merge(days,prcp_wide,by="CHARDATE",all=T)

##Computes a median across all stations for each day.
prcp_median <- apply(prcp_wide_all[,4:dim(prcp_wide_all)[2]],MARGIN=1,FUN=median,na.rm=T)

# ##Plots that along with the data.
# plot(log(prcp_dat$PRCP+0.01)~as.Date(prcp_dat$DATE),cex=0.1,col=rgb(0.5,0.5,0.5,0.5))
# lines(as.Date(days$date),log(prcp_median+0.01),lty=1)

##Eliminates all values more than 100mm from the median.
for (i in 4:dim(prcp_wide_all)[2]){
	diff <- abs(prcp_wide_all[,i] - prcp_median)
	out <- which(diff>1000)
	prcp_wide_all[out,i] <- NA
}

##Binds everything together and writes to disk.
colnames(tmax_wide_all) <- paste("TMAX",colnames(tmax_wide_all),sep="_")
colnames(tmin_wide_all) <- paste("TMIN",colnames(tmin_wide_all),sep="_")
colnames(prcp_wide_all) <- paste("PRCP",colnames(prcp_wide_all),sep="_")

alldat_wide <- data.frame(cbind(tmax_wide_all,tmin_wide_all,prcp_wide_all))
write.csv(alldat_wide,"stationdat_cleaned_1902_2013.csv")
