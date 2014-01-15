############################################## 

# TITLE: Snow cover algorithm


# CONTACT: Kevin Ford (krford@uw.edu) & Steve Kroiss (skroiss@gmail.com)


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

system.time({
    # start the clock to record how long the code takes to run
    
    
    ## READ IN THE DATA FILES FOR EACH TEMPERATURE SENSOR
    
    # For PC
    file.directory <- "C:/Users/Steve/Dropbox/Seeds&Seedlings(Steve)/R code - Seeds and Seedlings/Data - HOBO Spring 2013"
    setwd(file.directory)
    files <- list.files()  # the data files for each temperature sensor to be analyzed
    # calib.file <- read.csv('C:/Users/Steve/Dropbox/MORA Data & Logistics/MORA
    # CleanData/Microclimate_working/HoboIbuttonDatabase.csv', na.strings=c('', 'n/a')) # Import calibration file.  Replace
    # blank and n/a values with NA's
    
    figure.directory <- "C:/Users/Steve/Dropbox/Seeds&Seedlings(Steve)/R code - Seeds and Seedlings/Figures - Soil temp and snow cover"
    output.directory <- "C:/Users/Steve/Dropbox/Seeds&Seedlings(Steve)/R code - Seeds and Seedlings"
    
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
    
    for (k in 1:length(files)) {
        
        setwd(file.directory)
        d <- read.csv(files[k])
        d <- d[complete.cases(d[, 1:6]), ]
        
        # Name the columns
        names(d)[1] <- "YEAR"
        names(d)[2] <- "MONTH"
        names(d)[3] <- "DAY"
        names(d)[4] <- "HOUR"
        names(d)[5] <- "TEMP"
        
        ######################################## 
        
        # EXTRACT STAND INFO FROM FILENAME
        
        stand[k] <- strsplit(files[k], "-")[[1]][1]
        plot[k] <- strsplit(strsplit(files[k], "-S")[[1]][2], split = ".csv")[[1]][1]
        year[k] <- max(d$YEAR)
        
        
        ############################################## READING IN DATA FROM ONE FILE
        
        Temp <- as.numeric(d$TEMP)  #pulls out temperature data for each data row 
        
        # Pull out the date for each data row:
        d$Date <- as.Date(paste(d$MONTH, "/", d$DAY, "/", d$YEAR, sep = ""), format = "%m/%d/%Y")  #pulls out the date data for each row
        d$DOY <- as.numeric(format(d$Date, format = "%j"))  # find the unique days
        
        
        ############################################### 
        
        
        ## CRITERION 1
        RangeThresh <- 1  #the threshold temperature range in degrees C (i.e. if the range of soil temperatures on a given day exceed RangeThresh, snow probably absent)
        
        # Calculate the daily temperature range for each day:
        DOY <- d$DOY
        daily_range <- ave(d$TEMP, d$DOY, FUN = function(x) {
            diff(range(x))
        })
        daily_mean <- ave(d$TEMP, d$DOY, FUN = function(x) {
            (mean(x))
        })
        
        temp.d <- data.frame(d$Date, DOY, daily_range, daily_mean)
        d.unique <- unique(temp.d)
        head(d.unique)
        
        
        
        # Convert the daily temperature range into a binary code: 1 = daily temperature range did not exceed threshold (RangeThresh)
        # 0 = daily temperature range did exceed threshold (RangeThresh)
        
        d.unique$range_logical <- c()  # will store the binary code for daily temperature range for each day
        for (i in 1:length(d.unique$daily_range)) {
            d.unique$range_logical[i] <- ifelse(d.unique$daily_range[i] <= RangeThresh, yes = 1, no = 0)
        }
        
        
        ############################################## 
        
        
        # DETERMINE CALIBRATION TEMPERATURE
        calibration.temp = mean(d.unique$daily_mean[d.unique$range_logical == 1])  # calculate the mean temp for all the days that the temp didn't exceed RangeThresh
        calibration[k] <- ifelse(!is.na(calibration.temp), yes = calibration.temp, no = 0)  # If NA, set calibration = 0
        d$TEMP.calib <- d$TEMP - calibration[k]  # recalibrate data 
        
        
        ############################################### 
        
        
        ## CRITERION 2
        
        MaxThresh <- 2  # the threshold maximum temperature (i.e. if the max temperature on a given day exceeded MaxThresh, snow probably absent)
        
        # Calculate maximum daily temperature for each date
        maxT <- ave(d$TEMP.calib, d$DOY, FUN = function(x) {
            max(x)
        })
        
        max.d <- data.frame(DOY, maxT)
        max.unique <- unique(max.d)
        head(max.unique)
        
        d.unique$maxT <- max.unique$maxT
        
        
        # Convert the maximum daily temperature value into a binary code 1=daily max temperature did not exceed threshold
        # (MaxThresh) 0=daily max temperature did exceed threshold (MaxThresh)
        
        d.unique$maxT_logical <- c()  # will store the binary code for daily temperature range for each day
        for (i in 1:length(d.unique$maxT)) {
            d.unique$maxT_logical[i] <- ifelse(d.unique$maxT[i] <= MaxThresh, yes = 1, no = 0)
        }
        
        
        ################################################ 
        
        
        ## EVALUATE WHETHER OR NOT SNOW COVERED THE SENSOR ON EACH DATE BASED ON THE ABOVE TWO CRITERIA
        
        d.unique$snow_cover = c()  # will store the algorithm's evaluation of snow cover for each date (1=snow present, 0=snow absent)
        for (i in 1:length(d.unique$DOY)) {
            d.unique$snow_cover[i] <- ifelse(d.unique$range_logical[i] == 1 & d.unique$maxT_logical[i] == 1, yes = 1, no = 0)
        }
        
        d.unique <- d.unique[2:(length(d.unique$DOY) - 1), ]  #removes the first and last days, so there aren't any incomplete dates
        
        
        ################################################ 
        
        
        ## SUMMARIZING
        snow_appearance_date[k] <- as.character(d.unique$d.Date[1])  # first day when snow covered sensor
        snow_disappearance_date[k] <- as.character(d.unique$d.Date[dim(d.unique)[1]])  # last day when snow covered sensor
        snow_cover_duration[k] <- sum(d.unique$snow_cover)  #'snow cover duration' the total number of days with snow cover
        
        ################################################ 
        
        
        ## PLOT SOIL TEMPERATURE AND THE SNOW COVER ALGORITH OUTPUT TO MAKE SURE OUTPUT IS REASONABLE
        
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
        plot(d.unique$d.Date, d.unique$snow_cover, pch = 20, col = "blue", ylim = c(0, 1.2), axes = FALSE, xaxt = "n", yaxt = "n", 
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
    
    # Consolidate summarized results for each sensor into one data frame
    x <- 1:length(files)
    output <- data.frame(files[x], year[x], stand[x], plot[x], calibration[x], snow_appearance_date[x], snow_disappearance_date[x], 
        snow_cover_duration[x])
    
})  #stop the clock




# Save output file
setwd(output.directory)
write.table(output, file = "Snow_cover_summary_13.csv", sep = ",", row.names = FALSE)
 
