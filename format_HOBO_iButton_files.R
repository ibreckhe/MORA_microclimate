########## TEMPERATURE SENSOR DATA PROCESSING ##########

## by kevin ford (krford@uw.edu)

## The goal of this script is to take .csv files outputted by ibutton or HOBO temperature sensors and convert them to .csv
## files that can be inputted to R.

## Values in the '#===#' box below need to be changed for each run.  Once these values are changed, select all text and run.


# ============================================================================#

## give the year the data were collected, if necessary (will be appended to all processed file names)
year = "_2012-13"  # if year is not a part of file name, run this line of code (but not the line below), giving the appropriate year designation, to append year to each file name
# year = '' # if year is a part of file name, run this line of code (but not the line above) - this option will ensure that
# the year is not appended to processed file names

## name the working directories

# name the working directory with the .csv files you want to process.  there shouldn't be any .csv files in the below
# working directory that are not temperature sensor output files because the script will treat all .csv files as ibutton or
# HOBO files. There can be other types of files in the below working directory, this script will ignore them.
input.WorkingDirectory = "C:/Users/Kevin Ford/Documents/Science/projects/biome range shifts/analyses/seedling transplants/temperature_sensors/raw_2012-13"

# name the working directory where you want to put the processed .csv files
output.WorkingDirectory = "C:/Users/Kevin Ford/Documents/Science/projects/biome range shifts/analyses/seedling transplants/temperature_sensors/cleaned_2012-13"

# ============================================================================#



## set the working directory to the folder with all the raw .csv files you wish to process (this folder can contain other
## file types but should contain .csv files that are not HOBO output files)
setwd(input.WorkingDirectory)

## write a function for modifying the .csv files
modiFile <- function(CSV_FILE) {
    file <- read.csv(CSV_FILE, header = FALSE)  # read in the raw data .csv file as a dataframe
    
    sensor.test <- file[1, 1]
    is.ibutton <- grepl("iButton", sensor.test)  # determine whether the sensor is an ibutton or not
    
    if (is.ibutton == TRUE) 
        {
            # do this if the sensor is an ibutton do this if the ibutton is model DS1921G
            if (grepl("DS1921G", sensor.test) == TRUE) {
                file <- file[-c(1:15), ]  # remove rows from the top of the data frame
            }
            
            if (grepl("DS1922L", sensor.test) == TRUE) {
                # do this if the ibutton is model DS1922L
                file <- file[-c(1:20), ]  # remove rows from the top of the data frame
            }
            
            if (class(file) == "factor") {
                if (file[1] == "Date/Time") {
                  file = file[-c(1:3)]
                }
                
                if (file[1] == "Unit") {
                  file = file[-c(1:2)]
                }
                
                if (file[1] == "Value") {
                  file = file[-1]
                }
                
                file = as.character(file)
                tmp.file <- matrix(data = file, ncol = 3, nrow = length(file)/3, byrow = TRUE)  # take data (which is currently in vector format) and convert to matrix format
                file <- as.data.frame(tmp.file)  # convert to data frame
                file <- file[, -2]  # remove the second column (unit)
            }
            
            if (dim(file)[2] == 1) {
                file = as.character(file)
                tmp.file <- matrix(data = file, ncol = 3, nrow = length(file)/3, byrow = TRUE)  # take data (which is currently in vector format) and convert to matrix format
                file <- as.data.frame(tmp.file)  # convert to data frame
                file <- file[, -2]  # remove the second column (unit)
            }
            
            if (dim(file)[2] == 2) {
                file = file
            }
            
            if (dim(file)[2] == 3) {
                file <- file[, -2]  # remove the second column (unit)
            }
            
            light.fill <- array("", dim = c(dim(file)[1], 1))  # make a vector of NAs to fill the light column (ibuttons don't record light)
            file[, 3] = light.fill  # add the vector to the data frame
            
        }  # end ibutton-specific procedure
 else {
        # do this if the sensor is NOT an ibutton, and is a HOBO
        
        if (dim(file)[2] == 1) {
            # do this if the .csv file is tab-delimited (if tab-delimited, there will only be one column)
            file <- read.csv(CSV_FILE, header = F, sep = "\t")  # read in the raw data .csv file (as tab delimited) as a dataframe
        }
        
        file <- file[-c(1:2), ]  # remove the first two rows from the data frame
        file <- file[, 2:4]  # remove all columns except 2 and 3 (date/time and temperature)
    }  # end HOBO-specific procedure
    
    names(file) <- c("DateTime", "Temperature", "Light")  # name the columns
    return(file)
}

files <- list.files()  # lists all the files in the working directory
files <- files[which(grepl(".csv", files) == TRUE)]  # lists all the .csv files in the working directory (!!! WARNING: DON'T INCLUDE '.csv' IN ANY FILE NAMES IN THE WORKING DIRECTORY, '.csv' SHOULD ONLY APPEAR AS A FILE EXTENSION !!!). This script thinks that any file with the string '.csv' in the file name is a .csv file
dfs <- lapply(X = files, FUN = modiFile)  # modify all .csv files in the working directory


## append year to file names

locns = strsplit(files, ".csv")

appyear <- function(LOCNS) {
    l <- LOCNS
    sn <- as.character(c(l, year, ".csv"))
    newlocn <- paste(sn, collapse = "")
    return(newlocn)
}

locations <- sapply(X = locns, FUN = appyear)


## write the modified files as .csv files to an output folder

# set the working directory to the folder that will collect the output files (should be empty prior to running this script)
setwd(output.WorkingDirectory)

for (i in 1:length(locations)) {
    x = dfs[[i]]
    f = locations[i]
    write.csv(x, file = f, row.names = FALSE)
}

#  
