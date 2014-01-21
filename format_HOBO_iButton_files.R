########## TEMPERATURE SENSOR DATA PROCESSING ##########

## by kevin ford (krford@uw.edu)
## modified by Ian Breckheimer(ibreckhe@uw.edu)

## The goal of this script is to take .csv files outputted by ibutton or HOBO temperature sensors and convert them to .csv
## files that can be inputted to R.

## Values in the '#===#' box below need to be changed for each run.  Once these values are changed, select all text and run.


# ============================================================================#

## give the year the data were collected, if necessary (will be appended to all processed file names)
years <- c("-2012-13","-2011-12","-2011-12","-2010-11","-2010-11","-2009-10")  # if year is not a part of file name, run this line of code (but not the line below), giving the appropriate year designation, to append year to each file name
# year = '' # if year is a part of file name, run this line of code (but not the line above) - this option will ensure that
# the year is not appended to processed file names

## name the working directories

# name the working directory with the .csv files you want to process.  there shouldn't be any .csv files in the below
# working directory that are not temperature sensor output files because the script will treat all .csv files as ibutton or
# HOBO files. There can be other types of files in the below working directory, this script will ignore them.
input.WorkingDirectories <- c("~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2012-13",
                             "~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2011-12/ibuttondata2012",
                             "~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2011-12/HOBOdata2012",
                             "~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2010-11/dormant season",
                             "~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2010-11/growing season",
                             "~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2009-10/raw data")

# name the working directory where you want to put the processed .csv files
output.WorkingDirectory <-  "~/Dropbox/EcoForecasting_SDD_Phenology/Data&Analysis/Microclimate/compiled"

# prefix for output files. This should be a string that identifies the source data. e.g. ("Ford_HOBO_ibutton")
output.Prefixes <- c(rep("FordTheo",6))
# ============================================================================#



## set the working directory to the folder with all the raw .csv files you wish to process (this folder can contain other
## file types but should contain .csv files that are not HOBO output files)
setwd(input.WorkingDirectories[1])

## write a function for modifying the .csv files
modiFile <- function(CSV_FILE) {
    print(paste("Now processing file ",CSV_FILE))
    header <- scan(CSV_FILE,nlines=30,what='raw')  # read in the beginning of the file
    
    ## Tests to see what format the file is in.
    is.formatted <- header[1] %in% c("year,month,day,hour,temp,light",
                                     "year,month,day,hour,temp",
                                     "YEAR,MONTH,DAY,HOUR,TEMP,LIGHT",
                                     "YEAR,MONTH,DAY,HOUR,TEMP")
    is.datetime <- unlist(strsplit(header[1],split=","))[1] == "Date/Time"
    is.ibutton <- any(grepl("iButton", header))  # determine whether the sensor is an ibutton or not
    is.hobo <- any(grepl("GMT", header)) && any(grepl("File", header)) && any(grepl("Time", header))
    if (is.formatted == TRUE){
      file <- read.table(CSV_FILE,skip=1,sep=",")
    }
    else if (is.datetime) {
      # do this if the ibutton is model DS1922L
      file <- read.table(CSV_FILE,skip=1,header=FALSE,sep=",")  # remove rows from the top of the data frame
    }
    else if (is.ibutton == TRUE) {   
        # do this if the sensor is an ibutton do this if the ibutton is model DS1921G
            if (any(grepl("DS1921G", header)) == TRUE) {
                file <- read.table(CSV_FILE,skip=15,header=FALSE,sep=",")  # remove rows from the top of the data frame
            }
            if (any(grepl("DS1922L",header)) == TRUE) {
                # do this if the ibutton is model DS1922L
                file <- read.table(CSV_FILE,skip=20,header=FALSE,sep=",")  # remove rows from the top of the data frame
            }
            if (any(grepl("DS2422",header)) == TRUE) {
              # do this if the ibutton is model DS1922L
              file <- read.table(CSV_FILE,skip=20,header=FALSE,sep=",")  # remove rows from the top of the data frame
            }
            # do this if there is no header.
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
 else if (is.hobo==TRUE) {
      # do this if the sensor is NOT an ibutton, and is a HOBO
      sample <- header[c((length(header)-10):length(header))]
      sample_split <- unlist(strsplit(sample,split=","))
      is.comma.del <- length(sample)!=length(sample_split)
      if(is.comma.del) {
        file <- read.table(CSV_FILE, header = F, sep = ",",skip=2)  # read in the raw data .csv file (as comma delimited) as a dataframe
      }else{
        file <- read.table(CSV_FILE,header = F, sep="\t",skip=2)    # read in the raw data .csv file (as tab delimited) as a dataframe
      }
      file <- file[,2:4]  # remove all columns except 2 and 3 (date/time and temperature)
    }  # end HOBO-specific procedure
    names(file) <- c("DateTime", "Temperature", "Light")  # name the columns

    # Converts date vector to separate columns for year,month,day,and hour.
    dateTime <- strptime(file$DateTime, "%m/%d/%Y %H:%M")
    YEAR <- as.numeric(strftime(dateTime,format="%Y"))
    MONTH <- as.numeric(strftime(dateTime,format="%m"))
    DAY <- as.numeric(strftime(dateTime,format="%d"))
    HOUR <- as.numeric(strftime(dateTime,format="%H"))
    TEMP <- as.numeric(file$Temperature)
    file <- cbind(YEAR,MONTH,DAY,HOUR,TEMP)
    
    return(file)
}

for (i in 1:length(input.WorkingDirectories)){
  setwd(input.WorkingDirectories[i])
  print(paste("Now processing files in directory ",input.WorkingDirectories[i]))
  files <- list.files(".",pattern=".csv$") # lists all the .csv files in the working directory (!!! WARNING: DON'T INCLUDE '.csv' IN ANY FILE NAMES IN THE WORKING DIRECTORY, '.csv' SHOULD ONLY APPEAR AS A FILE EXTENSION !!!). This script thinks that any file with the string '.csv' in the file name is a .csv file
  dfs <- lapply(X = files, FUN = modiFile)  # modify all .csv files in the working directory
  
  
  ## append year to file names
  
  locns = strsplit(files, ".csv")
  
  appyear <- function(LOCNS) {
      l <- LOCNS
      sn <- as.character(c(l, years[i], ".csv"))
      newlocn <- paste(sn, collapse = "")
      return(newlocn)
  }
  
  locations <- sapply(X = locns, FUN = appyear)
  
  
  ## write the modified files as .csv files to an output folder
  
  # set the working directory to the folder that will collect the output files (should be empty prior to running this script)
  setwd(output.WorkingDirectory)
  
  for (j in 1:length(locations)) {
      x = dfs[[j]]
      f = paste(output.Prefixes[i],locations[j],sep="-")
      write.csv(x, file = f, row.names = FALSE)
  }
}

#  
