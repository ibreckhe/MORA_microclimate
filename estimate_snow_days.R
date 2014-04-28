############################################## 

# TITLE: Snow cover algorithm


# CONTACT: Kevin Ford (krford@uw.edu) & Steve Kroiss (skroiss@gmail.com) & Ian Breckheimer (ibreckhe@uw.edu)


# DESCRIPTION: The purpose of this code is to create an algorithm that determines whether snow is covering the ground on a
# particular date or not, based on temperature data from sensors at or just below the soil surface soil surface.  The
# algorithm was conceived and written in MATLAB by Mark Raleigh and then modified and translated into R by Kevin Ford and
# later modified by Steve Kroiss.

# Raleigh MS, Rittger K, Moore C, E., Henn B, Lutz J, A., et al. (2013) Ground-based testing of MODIS fractional snow cover
# in subalpine meadows and forests of the Sierra Nevada. Remote Sensing of Environment 128: 44-57.

# For the algorithm to consider a soil temperature sensor covered by snow, two criteria must be met: 1) the diurnal ground
# temperature range on that day does not exceed a certain value (Recommended: 1 deg C).

# 2) the temperature at that time step does not exceed a maximum value (Recommended: in theory, 0 deg C, but some records
# flatline as high as 2 deg C, likely due to bias in the sensor).

# In Mark's original algorithm, there was a third criterion: the temperature is constant for a user- specified number of
# timesteps. For Mark's original algorithm to score a sensor as 'snow covered' on a given date, two of the three criteria
# must be met. In Kevin's modified algorithm (this script), the algorithm only considers the first two criteria and both of
# those criteria must be met for the algorithm to score a sensor as 'snow covered' on a given date. Kevin did not use the
# third criterion because his dataset contained temperature sensors with different data logging intervals.

# You can try the recommended vaules for the algorithm criteria, but you should always check the output to make sure these
# parameters were appropriate. The idiosyncrasies of any one study require the algorithm to be modified to each study. This
# script plots soil temperature and snow algorithm output on the same graph for each sensor to help you do these checks.

# This script is an example of how to use the snow cover algorithm to estimate snow cover on a daily basis for two example
# temperature sensors (ibuttons deployed at Mount Rainier National Park).  You can modify the script to analyze your own
# data. However, the .csv data files containing your soil temperature measurements must have the exact same format as .csv
# files I provide as an example.


# NOTE ON CALIBRATION: The original script required a seperate file containing calibration data for each plot In this
# version, I have the code automatically determine a calibration value by averaging the terperature across all the days that
# meet Criteria 1 (<1 degree in temperature variation).  This appears to work very smoothly, but make sure to check the
# plots to verify that it worked correctly.


############################################## 

## READ IN THE DATA FILES FOR EACH TEMPERATURE SENSOR

##Loads required packages
library(data.table)

file.directory <- "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/compiled"
setwd(file.directory)
files <- list.files(pattern=".csv$")  # the data files for each temperature sensor to be analyzed
meta <- read.table("metadata.txt",sep=",",header=TRUE)

# calib.file <- read.csv('C:/Users/Steve/Dropbox/MORA Data & Logistics/MORA
# CleanData/Microclimate_working/HoboIbuttonDatabase.csv', na.strings=c('', 'n/a')) # Import calibration file.  Replace
# blank and n/a values with NA's

figure.directory <- "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/figs"
output.directory <- "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/processed/"

# # For MAC file.directory <- '~/Dropbox/Seeds&Seedlings(Steve)/R code - Seeds and Seedlings/Data - HOBO Spring 2013'
# setwd(file.directory) files <- list.files() # the data files for each temperature sensor to be analyzed # calib.file <-
# read.csv('C:/Users/Steve/Dropbox/MORA Data & Logistics/MORA CleanData/Microclimate_working/HoboIbuttonDatabase.csv', #
# na.strings=c('', 'n/a')) # Import calibration file.  Replace blank and n/a values with NA's # figure.directory <-
# 'C:/Users/Steve/Dropbox/Seeds&Seedlings(Steve)/R code - Seeds and Seedlings/Figures - Soil temp and snow cover' #
# output.directory <- 'C:/Users/Steve/Dropbox/Seeds&Seedlings(Steve)/R code - Seeds and Seedlings'

calibration <- c()
calibration.type <- c()
stand <- c()
plot <- c()
year <- c()
snow_appearance_date <- c()
snow_disappearance_date <- c()
snow_cover_duration <- c()  # in days

# start the clock to record how long the code takes to run
start_t <- Sys.time()

nfiles <- length(files)
for (k in 1:nfiles) {
    
    ##Prints progress.
    flush.console()
    print(paste("Now Processing file: ",files[k],"(",k," of",nfiles,")"))
  
    ############################################## READING IN DATA FROM ONE FILE
    
    setwd(file.directory)
    d <- read.csv(files[k])
    d <- d[complete.cases(d[, 1:5]), ]
    
    # Name the columns
    names(d)[1] <- "YEAR"
    names(d)[2] <- "MONTH"
    names(d)[3] <- "DAY"
    names(d)[4] <- "HOUR"
    names(d)[5] <- "TEMP"
    
    d$Date <- as.Date(paste(d$MONTH, "/", d$DAY, "/", d$YEAR, sep = ""), format = "%m/%d/%Y")
    d$DOY <- as.numeric(format(d$Date, format = "%j"))  # find the unique days
    d <- data.table(d,key="Date")
    
    ######################################## 
    
    # EXTRACT STAND INFO FROM FILENAME
    
    stand[k] <- strsplit(files[k], "_")[[1]][2]
    plot[k] <- strsplit(files[k], "_")[[1]][3]
    year[k] <- max(d$YEAR)
        
    ############################################### 
    
    ## Criteria
    RangeThresh <- 1  #the threshold temperature range in degrees C (i.e. if the range of soil temperatures on a given day exceed RangeThresh, snow probably absent)
    MaxThresh <- 2
    
    # Create an empty data table to store values
    days <- unique(d$Date)
    ndays <- length(days)
    
    daily <- data.table(date=unique(d$Date),
                        range=rep(NA,ndays),
                        mean=rep(NA,ndays),
                        rangethresh=rep(NA,ndays),
                        maxthresh=rep(NA,ndays),
                        snow=rep(NA,ndays))
    setkey(daily,"date")
                        
    # Calculate the mean daily temperature and temperature range for each day:
    daily[[2]] <- d[,diff(range(TEMP)),by=Date][[2]]
    daily[[3]] <- d[,mean(TEMP,na.rm=T),by=Date][[2]]
    daily[[4]] <- d[,diff(range(TEMP)) < RangeThresh,by=Date][[2]]
    daily[[5]] <- d[,max(TEMP) < MaxThresh,by=Date][[2]]

    ############################################## 
    
    # DETERMINE CALIBRATION TEMPERATURE
    calibration.temp <- mean(daily$mean[daily$rangethresh == TRUE])  # calculate the mean temp for all the days that the temp didn't exceed RangeThresh
    calibration[k] <- ifelse(!is.na(calibration.temp), yes = calibration.temp, no = 0)  # If NA, set calibration = 0
    d$TEMP.calib <- d$TEMP - calibration[k]  # recalibrate data 
    
    ################################################ 
  
    ## EVALUATE WHETHER OR NOT SNOW COVERED THE SENSOR ON EACH DATE BASED ON THE ABOVE TWO CRITERIA
    
    daily[[6]] <- daily[[4]] & daily[[5]]  # will store the algorithm's evaluation of snow cover for each date (1=snow present, 0=snow absent)
    d.snow <- subset(daily,snow==TRUE)
    ################################################ 
    
    ## SUMMARIZING
    snow_appearance_date[k] <- as.character(min(d.snow$date))  # first day when snow covered sensor
    snow_disappearance_date[k] <- as.character(max(d.snow$date))  # last day when snow covered sensor
    snow_cover_duration[k] <- sum(d.snow$snow)  #'snow cover duration' the total number of days with snow cover
  
    
    ## PLOT SOIL TEMPERATURE AND THE SNOW COVER ALGORITHm OUTPUT TO MAKE SURE OUTPUT IS REASONABLE
    
    # Save figure as pdf
    setwd(figure.directory)
    population <- strsplit(files[k], ".csv")[[1]]
    graph.file <- paste(population, ".pdf", sep = "")
    pdf(file = graph.file, width = 10, height = 7)
    
    # Left Y axis
    par(mar = c(4, 6, 3, 6))
    plot(d$Date, d$TEMP.calib, main = population, axes = F, xlab = "", ylab = "", pch = 20, col = "red")
    axis(2, col = "red", col.axis = "red", col.ticks = "red")
    leftY <- expression(paste("Temperature", degree, "C"))
    text(par("usr")[1] - 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj = 0.5, leftY, srt = 90, xpd = TRUE, 
        col = "red")
    
    # Right Y axis
    par(new = TRUE)
    plot(daily$date, daily$snow, pch = 20, col = "blue", ylim = c(0, 1.2), axes = FALSE, xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "")
    axis(4, col = "blue", col.axis = "blue", col.ticks = "blue", yaxp = c(0, 1, 1), labels = F)
    rightY <- paste("Snow cover")
    text(par("usr")[2] + 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj = 0.5, rightY, srt = 270, xpd = TRUE, 
        col = "blue")
    text(par("usr")[2] + 15, 0, 0, srt = 270, xpd = T, col = "blue")
    text(par("usr")[2] + 15, 1, 1, srt = 270, xpd = T, col = "blue")
    
    # X axis
    axis.Date(1, d$Date, at = seq(from = min(d$Date), to = max(d$Date), by = "months"))
    
    dev.off()
    
    ################################################ 
    
}

##Computes elapsed time.
stop_t <- Sys.time()
print(paste("Run took ",stop_t - start_t))

# Consolidate summarized results for each sensor into one data frame
output <- data.frame(files, 
                     year, 
                     stand, 
                     plot, 
                     calibration, 
                     snow_appearance_date, 
                     snow_disappearance_date,
                     snow_cover_duration)

##Match the snow cover information to the sensor metadata.
out_merged <- merge(meta,output,by.x="out_filename",by.y="files",all.x=T)

##Data quality flags:
out_merged$flag_sensor_fail <- (as.Date(out_merged$snow_disappearance_date) - as.Date(out_merged$date_max)) >= -1
out_merged$flag_temp_high <- out_merged$temp_max > 90
out_merged$flag_temp_low <- out_merged$temp_min < -30
out_merged$flag_high_calib <- out_merged$calibration >= 2 | out_merged$calibration <= -2
out_merged$flag_no_snow <- out_merged$snow_cover_duration <= 15
out_merged$flagged <- with(out_merged, flag_sensor_fail | flag_temp_high | flag_temp_low | flag_high_calib | flag_no_snow)

##Moves .csv and .pdf graphics of flagged files to a new directories
out_flagged <- out_merged[out_merged$flagged==TRUE,]
flagged_csvs <- out_flagged$out_filename
setwd(output.directory)
dir.create("./flagged")
file.copy(paste(file.directory,"/",as.character(flagged_csvs),sep=""),"./flagged")
flagged_pdfs <- sub(".csv",".pdf",out_flagged$out_filename)
setwd(figure.directory)
dir.create("./flagged")
file.copy(flagged_pdfs,"./flagged")

##Moves .csv and .pdf graphics of unflagged files to a new directory
out_unflagged <- out_merged[out_merged$flagged==FALSE,]
unflagged_csvs <- out_unflagged$out_filename
setwd(output.directory)
dir.create("./unflagged")
file.copy(paste(file.directory,"/",as.character(unflagged_csvs),sep=""),"./unflagged")
unflagged_pdfs <- sub(".csv",".pdf",out_unflagged$out_filename)
setwd(figure.directory)
dir.create("./unflagged")
file.copy(unflagged_pdfs,"./unflagged")

##subsets data for the unflagged files
out_unflagged<- out_merged[out_merged$flagged==FALSE,]

# Save output file
setwd(output.directory)
write.table(out_unflagged, file = "Snow_cover_meta_cleaned.csv", sep = ",", row.names = FALSE)
 
