## Script to play around with downloading climate data from the world bank.
library(devtools)
library(rnoaa)
library(plyr)
library(reshape2)

## Sets the API key for downloading NOAA data.
options(noaakey = "NjuEQnGAERfxLrZNIiHqMUJsSeqiPuRR")

## Sets the working directory.
setwd("~/code/mora_climate")

## Current stations in the vicinity of Mt. Rainier.  stationIDs <- c('GHCND:USS0021C38S', 'GHCND:USS0021B63S',
## 'GHCND:USS0021B13S', 'GHCND:USR0000WGRN', 'GHCND:USS0021B62S', 'GHCND:USC00454764', 'GHCND:US1WALW001',
## 'GHCND:USS0021C17S', 'GHCND:USS0021C40S', 'GHCND:USC00455704', 'GHCND:USC00456262', 'GHCND:USS0021C35S',
## 'GHCND:USS0021C33S', 'GHCND:USC00456896', 'GHCND:USC00456898', 'GHCND:USC00456909', 'GHCND:USS0021C28S')

## All current stations within 50km of the summit.
stationIDs <- c("GHCND:USS0021C35S", "GHCND:USS0021B62S", "GHCND:USC00455704", "GHCND:USS0021B63S", "GHCND:USC00456896", "GHCND:USS0021C28S", 
    "GHCND:USC00455224", "GHCND:USC00456909", "GHCND:USS0021C33S", "GHCND:USC00456898", "GHCND:USS0021C40S", "GHCND:USC00456262", 
    "GHCND:USS0021B13S", "GHCND:USC00454764", "GHCND:USS0021C17S", "GHCND:USS0021C38S", "GHCND:US1WALW0015", "GHCND:USS0021B42S", 
    "GHCND:USR0000WGRN", "GHCND:US1WAPR0028")

# ##All (current and historical) stations in the vicinity of Mt. Rainier.  stationIDs <- c('GHCND:USC00450285',
# 'GHCND:USC00450945', 'GHCND:USC00450969', 'GHCND:USS0021C38S', 'GHCND:USS0021B63S', 'GHCND:USC00451113',
# 'GHCND:USC00451110', 'GHCND:USS0021B13S', 'GHCND:USC00452493', 'GHCND:USC00452722', 'GHCND:USC00453177',
# 'GHCND:USC00454286', 'GHCND:USC00453219', 'GHCND:USR0000WGRN', 'GHCND:USC00453357', 'GHCND:USR0000WHAG',
# 'GHCND:USS0021B62S', 'GHCND:USC00454620', 'GHCND:USC00454619', 'GHCND:USC00456894', 'GHCND:USC00454764',
# 'GHCND:US1WALW0015', 'GHCND:USC00455425', 'GHCND:USS0021C17S', 'GHCND:USS0021C40S', 'GHCND:USC00455704',
# 'GHCND:USC00456201', 'GHCND:USC00456265', 'GHCND:USC00456262', 'GHCND:USS0021C35S', 'GHCND:USC00456381',
# 'GHCND:USS0021C33S', 'GHCND:USC00456892', 'GHCND:USC00456896', 'GHCND:USC00456898', 'GHCND:USC00456900',
# 'GHCND:USC00456909', 'GHCND:USC00458653', 'GHCND:USS0021C28S', 'GHCND:USC00459170', 'GHCND:USC00456385',
# 'GHCND:USC00459171')

## Function to get all daily temperature and precipitation data for the stations.  NOAA limits API calls to 1000 records, so
## this is a workaround.
get_ghcn_daily <- function(start_date, end_date, type, stations, chunk_days = 20) {
    dates <- seq(start_date, end_date, by = chunk_days)
    print("Downloading data for date chunk:")
    data_list <- list()
    for (i in 2:length(dates)) {
        print(paste(i," of ",length(dates),sep=""))
        data_list[[i - 1]] <- ncdc(datasetid = "GHCND", datatypeid = type, stationid = stations, startdate = as.character(dates[i - 
            1]), enddate = as.character(dates[i]), limit = 1000)
    }
    data <- data_list[[1]]
    for (i in 2:length(data_list)) {
        data <- ncdc_combine(data, data_list[[i]])
    }
    return(data)
}

## Gets the most recent data.
start <- as.Date("1981-01-01")
end <- as.Date("2015-12-31")
tmax_data <- get_ghcn_daily(start, end, type = "TMAX", stations = stationIDs, chunk_days = 20)[[1]]
tmin_data <- get_ghcn_daily(start, end, type = "TMIN", stations = stationIDs, chunk_days = 20)[[1]]
prcp_data <- get_ghcn_daily(start, end, type = "PRCP", stations = stationIDs, chunk_days = 20)[[1]]

## Quick plot to check validity.
par(mfrow = c(3, 1))
scatter.smooth(as.Date(tmax_data$date), tmax_data$value, pch = 20, cex = 0.2, col = rgb(0.5, 0.5, 0.5, 0.5), main = "TMAX", 
    lpars = list(col = "red"), span = 0.1)
scatter.smooth(as.Date(tmin_data$date), tmin_data$value, pch = 20, cex = 0.2, col = rgb(0.5, 0.5, 0.5, 0.5), main = "TMIN", 
    lpars = list(col = "red"), span = 0.1)
scatter.smooth(as.Date(prcp_data$date), prcp_data$value, pch = 20, cex = 0.2, col = rgb(0.5, 0.5, 0.5, 0.5), main = "PRCP", 
    lpars = list(col = "red"), span = 0.05)


## Writes the data to disk.
write.csv(tmax_data, "MORA_daily_climate_tmax_2013_2015.csv")
write.csv(tmin_data, "MORA_daily_climate_tmin_2013_2015.csv")
write.csv(prcp_data, "MORA_daily_climate_prcp_2013_2015.csv")

tmax_data <- read.csv("MORA_daily_climate_tmax_2013_2015.csv")
tmin_data <- read.csv("MORA_daily_climate_tmin_2013_2015.csv")
prcp_data <- read.csv("MORA_daily_climate_prcp_2013_2015.csv")


#### Cleans the contemporary climate data.####
tmax_data$date <- strptime(tmax_data$date,format="%Y-%m-%dT%H:%M:%S")
tmin_data$date <- strptime(tmin_data$date,format="%Y-%m-%dT%H:%M:%S")
prcp_data$date <- strptime(prcp_data$date,format="%Y-%m-%dT%H:%M:%S")

##TMAX cleaning.
tmax_data <- tmax_data[,c("date","station","value")]
tmax_data <- unique(tmax_data)
tmax_data$value[tmax_data$value > 500] <- NA
tmax_data$value[tmax_data$value < -300] <- NA

## Creates a continuous sequence of days.
start <- min(tmax_data$date)
end <- max(tmax_data$date)
days <- data.frame(date = seq(start, end, by = "day"))

tmax_wide <- dcast(tmax_data, date ~ station, value.var = "value")
dim(tmax_wide)

## Eliminates all values more than 30C from the median.
tmax_median <- apply(tmax_wide[, 2:dim(tmax_wide)[2]], MARGIN = 1, FUN = median, na.rm = T)
for (i in 2:dim(tmax_wide)[2]) {
  diff <- abs(tmax_wide[, i] - tmax_median)
  out <- which(diff > 300)
  out2 <- which(diff < -300)
  tmax_wide[out, i] <- NA
  tmax_wide[out2, i] <- NA
}

##TMIN cleaning.
tmin_data <- tmin_data[,c("date","station","value")]
tmin_data <- unique(tmin_data)
tmin_data$value[tmin_data$value > 500] <- NA
tmin_data$value[tmin_data$value < -300] <- NA

## Creates a continuous sequence of days.
start <- min(tmin_data$date)
end <- max(tmin_data$date)
days <- data.frame(date = seq(start, end, by = "day"))

tmin_wide <- dcast(tmin_data, date ~ station, value.var = "value")
dim(tmin_wide)

## Eliminates all values more than 30C from the median.
tmin_median <- apply(tmin_wide[, 2:dim(tmin_wide)[2]], MARGIN = 1, FUN = median, na.rm = T)
for (i in 2:dim(tmin_wide)[2]) {
  diff <- abs(tmin_wide[, i] - tmin_median)
  out <- which(diff > 300)
  out2 <- which(diff < -300)
  tmin_wide[out, i] <- NA
  tmin_wide[out2, i] <- NA
}

##PRCP cleaning.
prcp_data <- prcp_data[,c("date","station","value")]
prcp_data <- unique(prcp_data)

## Creates a continuous sequence of days.
start <- min(prcp_data$date)
end <- max(prcp_data$date)
days <- data.frame(date = seq(start, end, by = "day"))

prcp_wide <- dcast(prcp_data, date ~ station, value.var = "value")
dim(prcp_wide)

## Eliminates all values more than 300mm from the median.
prcp_median <- apply(prcp_wide[, 2:dim(prcp_wide)[2]], MARGIN = 1, FUN = median, na.rm = T)
for (i in 2:dim(prcp_wide)[2]) {
  diff <- abs(prcp_wide[, i] - prcp_median)
  out <- which(diff > 300)
  out2 <- which(diff < -300)
  prcp_wide[out, i] <- NA
  prcp_wide[out2, i] <- NA
}

## Binds everything together and writes to disk.
colnames(tmax_wide) <- c("date",paste("TMAX", colnames(tmax_wide)[-1], sep = "_"))
colnames(tmin_wide) <- c("date",paste("TMIN", colnames(tmin_wide)[-1], sep = "_"))
colnames(prcp_wide) <- c("date",paste("PRCP", colnames(prcp_wide)[-1], sep = "_"))

recentdat_wide <- data.frame(cbind(tmax_wide, tmin_wide[,-1], prcp_wide[,-1]))
write.csv(recentdat_wide,"tmax_2013_2015_wide.csv",row.names=FALSE)


#### Brings in and cleans the historical station data.####
statdat <- read.csv("station_data_50km_1981_2013.csv")

## Cleans the TMAX data.
tmax_dat <- statdat[, c("STATION", "LATITUDE", "LONGITUDE", "DATE", "TMAX")]
tmax_dat$DATE <- strptime(as.character(tmax_dat$DATE), format = "%Y%m%d")

# ##Examines the data.  scatter.smooth(as.Date(tmax_dat$DATE),tmax_dat$TMAX,pch=20,cex=0.2,
# col=rgb(0.5,0.5,0.5,0.5),main='TMAX',lpars=list(col='red'))

## Eliminates impossible values
tmax_dat$TMAX[which(tmax_dat$TMAX > 500)] <- NA
tmax_dat$TMAX[which(tmax_dat$TMAX < -300)] <- NA

## Creates a continuous sequence of days.
start <- min(tmax_dat$DATE)
end <- max(tmax_dat$DATE)
days <- data.frame(date = seq(start, end, by = "day"))

## Casts the tmax data to 'wide' format.
tmax_wide <- dcast(tmax_dat, DATE ~ STATION, value.var = "TMAX")
dim(tmax_wide)

## Computes a median across all stations for each day.
tmax_median <- apply(tmax_wide[, 2:dim(tmax_wide)[2]], MARGIN = 1, FUN = median, na.rm = T)

## Eliminates all values more than 15C from the median.
for (i in 2:dim(tmax_wide)[2]) {
    diff <- abs(tmax_wide[, i] - tmax_median)
    out <- which(diff > 150)
    tmax_wide[out, i] <- NA
}

#### Cleans the TMIN data.####
tmin_dat <- statdat[, c("STATION", "LATITUDE", "LONGITUDE", "DATE", "TMIN")]
tmin_dat$DATE <- strptime(as.character(tmin_dat$DATE), format = "%Y%m%d")

## Examines the data.  scatter.smooth(as.Date(tmin_dat$DATE),tmin_dat$TMIN,pch=20,cex=0.2,
## col=rgb(0.5,0.5,0.5,0.5),main='TMIN',lpars=list(col='red'))

## Eliminates impossible values
tmin_dat$TMIN[which(tmin_dat$TMIN > 250)] <- NA
tmin_dat$TMIN[which(tmin_dat$TMIN < -300)] <- NA

## Creates a continuous sequence of days.
start <- min(tmin_dat$DATE)
end <- max(tmin_dat$DATE)
days <- data.frame(date = seq(start, end, by = "day"))

## Casts the tmax data to 'wide' format.
tmin_wide <- dcast(tmin_dat, DATE ~ STATION, value.var = "TMIN")
dim(tmin_wide)

## Computes a median across all stations for each day.
tmin_median <- apply(tmin_wide[, 2:dim(tmin_wide)[2]], MARGIN = 1, FUN = median, na.rm = T)

## Eliminates all values more than 15C from the median.
for (i in 2:dim(tmin_wide)[2]) {
    diff <- abs(tmin_wide[, i] - tmin_median)
    out <- which(diff > 150)
    tmin_wide[out, i] <- NA
}

#### Cleans the PRCP data.####
prcp_dat <- statdat[, c("STATION", "LATITUDE", "LONGITUDE", "DATE", "PRCP")]
prcp_dat$DATE <- strptime(as.character(prcp_dat$DATE), format = "%Y%m%d")

## Examines the data.
scatter.smooth(as.Date(prcp_dat$DATE), prcp_dat$PRCP, pch = 20, cex = 0.2, col = rgb(0.5, 0.5, 0.5, 0.5), main = "PRCP", lpars = list(col = "red"))

## Eliminates impossible values
prcp_dat$PRCP[which(prcp_dat$PRCP > 2000)] <- NA
prcp_dat$PRCP[which(prcp_dat$PRCP < 0)] <- NA

## Casts the prcp data to 'wide' format.
prcp_wide <- dcast(prcp_dat, DATE ~ STATION, value.var = "PRCP")
dim(prcp_wide)

## Computes a median across all stations for each day.
prcp_median <- apply(prcp_wide[, 2:dim(prcp_wide)[2]], MARGIN = 1, FUN = median, na.rm = T)

# ##Plots that along with the data.  plot(log(prcp_dat$PRCP+0.01)~as.Date(prcp_dat$DATE),cex=0.1,col=rgb(0.5,0.5,0.5,0.5))
# lines(as.Date(days$date),log(prcp_median+0.01),lty=1)

## Eliminates all values more than 100mm from the median.
for (i in 2:dim(prcp_wide)[2]) {
    diff <- abs(prcp_wide[, i] - prcp_median)
    out <- which(diff > 1000)
    prcp_wide[out, i] <- NA
}

## Binds everything together and writes to disk.
colnames(tmax_wide) <- c("date",paste("TMAX", colnames(tmax_wide)[-1], sep = "_"))
colnames(tmin_wide) <- c("date",paste("TMIN", colnames(tmin_wide)[-1], sep = "_"))
colnames(prcp_wide) <- c("date",paste("PRCP", colnames(prcp_wide)[-1], sep = "_"))

alldat_wide <- data.frame(cbind(tmax_wide, tmin_wide[,-1], prcp_wide[,-1]))
write.csv(alldat_wide, "stationdat_cleaned_1981_2013.csv")

##Joins recent with older data.
hist_recent_wide <- merge(alldat_wide,recentdat_wide,all=TRUE)
write.csv(hist_recent_wide,"GHCN_daily_1981_2015.csv",row.names=FALSE)

##Plots data.
start <- as.POSIXct("2001-10-01")
end <- as.POSIXct("2001-12-07")
plot(hist_recent_wide$date,hist_recent_wide$TMIN_GHCND.USC00454764,type="l",xlim=c(start,end))
points(hist_recent_wide$date,hist_recent_wide$TMIN_GHCND.USC00455224,type="l",col=2)
points(hist_recent_wide$date,hist_recent_wide$TMIN_GHCND.USC00455704,type="l",col=3)
points(hist_recent_wide$date,hist_recent_wide$TMIN_GHCND.USC00456262,type="l",col=4)
points(hist_recent_wide$date,hist_recent_wide$TMIN_GHCND.USC00456896,type="l",col=5)
points(hist_recent_wide$date,hist_recent_wide$TMIN_GHCND.USC00456898,type="l",col=6)
