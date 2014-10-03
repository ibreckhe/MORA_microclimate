########## TEMPERATURE SENSOR DATA PROCESSING ##########

## by kevin ford (krford@uw.edu)
## modified by Ian Breckheimer(ibreckhe@uw.edu)

## The goal of this script is to take .csv files outputted by ibutton or HOBO temperature sensors and convert them to .csv
## files that can be inputted to R.

## Values in the '#===#' box below need to be changed for each run.  Once these values are changed, select all text and run.


#### ==================INPUTS======================================####

## give the year the data were collected, if necessary (will be appended to all processed file names)
# if year is not a part of file name, run this line of code (but not the line below), giving the appropriate year designation, to append year to each file name
# year = '' # if year is a part of file name, run this line of code (but not the line above) - this option will ensure that
# the year is not appended to processed file names

## name the working directories

# name the working directories with the .csv files you want to process.
input.WorkingDirectories <- c("~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2012-13",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2011-12/ibuttondata2012",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2011-12/HOBOdata2012",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2010-11/dormant season",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2010-11/growing season",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Elli/temperature sensors/2009-10/manipulated",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/hobodata2013",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/ibuttondata2012",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/HoboData2012Summer",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/HoboData2012Fall",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/HoboData2011Summer",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/HoboData2011Fall",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/HoboData2010-06",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ailene/HoboIbuttonData/HoboData2009-09",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2013/air_iButton",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2013/air_HOBO/summerDownload",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2013/air_HOBO/springDownload",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2012",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2011-2012",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2010-2011",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2009-2010",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/AirTemp&Humidity/cleaned/2008-2009",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2013/SummerDownload_Clean",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2013/SpringDownload_Clean",                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2013/SummerDownload_Clean",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2012/Spring-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2012/Fall-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2012/Spring-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2011/Fall-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2011/Spring-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2010/Fall-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2010/Spring-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2009/Fall-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/SoilTemp&Light/Clean_Data/2009/Spring-HOBO",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/2014/Air Hobo (Summer Download)",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/2014/Humidity Ibutton (Summer Download)",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/2014/Soil Hobo (Summer Download)",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/2014/Wonderland 2014",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/Wonderland/2013",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/Wonderland/2012",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/Wonderland/2011",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/Wonderland/2010",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/JHRL/Wonderland/2009",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Jessica",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ian/Riparian 2014/spring",
                              "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/Ian/Riparian 2014/fall")

                        
# name the working directory where you want to put the processed .csv files
output.WorkingDirectory <-  "~/Dropbox/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/compiled_air"

# prefix for output files. This should be a string that identifies the source data. e.g. ("Ford_HOBO_ibutton")
output.Prefixes <- c(rep("FordTheo",6),
                     rep("Ettinger",8),
                     rep("JHRL",29),
                     rep("Lundquist",1),
                     rep("Breckheimer",2))

## set the working directory to the folder with all the raw .csv files you wish to process (this folder can contain other
## file types but should contain .csv files that are not HOBO output files)
setwd(input.WorkingDirectories[1])


#### ==================FUNCTIONS======================================####

## write a function for modifying the .csv files
formatMicro <- function(CSV_FILE) {
    print(paste("Now processing file ",CSV_FILE))
    header <- scan(CSV_FILE,nlines=30,what='raw',quiet=TRUE)  # read in the beginning of the file
    
    ## Tests to see what format the file is in.
    is.formatted <- tolower(header[1]) %in% c("year,month,day,hour,temp,light",
                                              "year,month,day,hour,temperature",
                                              "year,month,day,hour,temp",
                                              "year,month,day,hour,temp,",
                                              "year,month,day,hour,temp,,",
                                              "year,month,day,hour,temp,,,",
                                              "year,month,day,hour,temp,,,,",
                                              "year,month,day ,hour,temp",
                                              "year,month,day,hour,value",
                                              "year,month,day,hour,temp,light,,",
                                              "year,month,day,hour,temp,light,,,",
                                              "year,month,day,hour,temp,light,,,,",
                                              "year,month,day,hour,temp,light,,,,,,,",
                                              "year,month,day,hour,temp,intensity",
                                              "year,month,day,hour,temp,temp")
                                              
    is.datetime <- unlist(strsplit(header[1],split=","))[1] == "Date/Time"
    is.ibutton <- any(grepl("iButton", header))  # determine whether the sensor is an ibutton or not
    is.hobo <- any(grepl("File", header)) && any(grepl("Time", header)) && any(grepl("End", header))
    if (is.formatted == TRUE){
      file <- read.table(CSV_FILE,skip=1,sep=",")
      DateTime <- paste(file[,2],"/",file[,3],"/",file[,1]," ",file[,4],":00",sep="")
      file <- cbind(DateTime,file[-c(1:4)],NA)
    }
    else if (is.datetime) {
      table <- read.table(CSV_FILE,skip=1,header=FALSE,sep=",")  # remove rows from the top of the data frame
      DateTime <- strptime(table[,1],format="%m/%d/%y %H:%M")
      temp <- as.vector(table[,2])
      file <- data.frame(DateTime,temp,light=NA)
    }
    else if (is.ibutton == TRUE) {   
        # do this if the sensor is an ibutton do this if the ibutton is model DS1921G
            if (any(grepl("DS1921G", header)) == TRUE) {
                file <- read.table(CSV_FILE,skip=15,header=FALSE,sep=",")  # remove rows from the top of the data frame
            }
            if (any(grepl("DS1922L",header)) == TRUE) {
                # do this if the ibutton is model DS1922L
                file <- read.table(CSV_FILE,skip=22,header=FALSE,sep=",")  # remove rows from the top of the data frame
            }
            if (any(grepl("DS2422",header)) == TRUE) {
              # do this if the ibutton is model DS1922
              file <- read.table(CSV_FILE,skip=20,header=FALSE,sep=",")  # remove rows from the top of the data frame
            }
            if (any(grepl("DS1923",header)) == TRUE) {
              # do this if the ibutton is model DS1923
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
        file <- read.table(CSV_FILE, header = F, sep = ",",skip=2,fill=TRUE)  # read in the raw data .csv file (as comma delimited) as a dataframe
      }
      else{
        file <- read.table(CSV_FILE,header = F, sep="\t",skip=2,fill=TRUE)    # read in the raw data .csv file (as tab delimited) as a dataframe
      }
      file <- file[,2:4]  # remove all columns except 2 and 3 (date/time and temperature)
    }  # end HOBO-specific procedure
  else{
    print(paste("Could not recognize the formatting of ",CSV_FILE))
  }
  
  names(file) <- c("DateTime", "Temperature", "Light")  # name the columns
  
  # Attempts to extract the time-zone from the header.
  if(any(grepl("GMT-07:00", header,fixed=T)) | any(grepl("PDT", header,fixed=T))){
    tz <- "PDT"
  }
  else if(any(grepl("GMT-00:00", header,fixed=T))){
    tz <- "GMT"
  }
  else{
    tz <- "UNK"
  }
  
  # Tests for a valid 4-digit year
  valid_years <- 2007:2014
  datestring <- strsplit(as.character(file$DateTime[1]),split=" ")[[1]][1]
  yeartest <- function(x){any(grepl(as.character(x),datestring))}
  valid.year <- any(sapply(valid_years,FUN=yeartest))
  am.pm <- any(grepl("PM",file$DateTime[1:20]))
    
  # Converts date vector to separate columns for year,month,day,and hour.
  if(class(file$DateTime)[1]=="POSIXct"){
    dateTime <- file$DateTime
  }else if(am.pm && valid.year) {
      dateTime <- strptime(file$DateTime, "%m/%d/%Y %r")
  }else if(am.pm==FALSE && valid.year==TRUE) {
    dateTime <- strptime(file$DateTime, "%m/%d/%Y %H:%M")
  }else if(am.pm==TRUE && valid.year==FALSE) {
    dateTime <- strptime(file$DateTime, "%m/%d/%y %r")
  }else{
    dateTime <- strptime(file$DateTime, "%m/%d/%y %H:%M")
  }
    
  YEAR <- as.numeric(strftime(dateTime,format="%Y"))
  MONTH <- as.numeric(strftime(dateTime,format="%m"))
  DAY <- as.numeric(strftime(dateTime,format="%d"))
  HOUR <- as.numeric(strftime(dateTime,format="%H"))
  MIN <- as.numeric(strftime(dateTime,format="%M"))
  TEMP <- as.numeric(file$Temperature)
  file <- cbind(YEAR,MONTH,DAY,HOUR,MIN,TEMP)
  
  # Extracts the range of dates in the file.
  date_min <- min(dateTime,na.rm=T)
  date_max <- max(dateTime,na.rm=T)
  date_lab <- paste(strftime(date_min,format="%Y-%m"),strftime(date_max,format="%Y-%m"),sep="-")

  # Measures the logging interval
  log_int <- dateTime[2] - dateTime[1]
    
  # Extracts the minimum and maximum temperature.
  temp_max <- max(TEMP,na.rm=T)
  temp_min <- min(TEMP,na.rm=T)
    
  # Extracts the number of measurements.
  n_measurements <- length(na.omit(TEMP))

  file_attrib <- list(data=file,
                      filename=strsplit(CSV_FILE,split=".csv")[[1]],
                      filepath=paste(getwd(),CSV_FILE,sep="/"),
                      n_measurements=n_measurements,
                      log_interval=log_int,
                      date_min=date_min,
                      date_max=date_max,
                      date_lab=date_lab,
                      temp_min=temp_min,
                      temp_max=temp_max,
                      hobo=is.hobo,
                      ibutton=is.ibutton,
                      formatted=is.formatted,
                      datetime=is.datetime,
                      tz=tz)
  
  return(file_attrib)
}


#### ==================FILE PROCESSING==============================####

##Checks to see how many input files there are.
csv_all <- c()
for (i in 1:length(input.WorkingDirectories)){
  csv_files <- list.files(input.WorkingDirectories[i],pattern=".csv$")
  csv_all <- c(csv_all,csv_files)
}
nfiles <- length(csv_all)
print(paste("Now processing ",nfiles," microclimate files."))

##Sets up a file counter.
file_n <- 0

for (i in 1:length(input.WorkingDirectories)){
  setwd(input.WorkingDirectories[i])
  print(paste("Now processing files in directory ",input.WorkingDirectories[i]))
  files <- list.files(".",pattern=".csv$") # lists all the .csv files in the working directory (!!! WARNING: DON'T INCLUDE '.csv' IN ANY FILE NAMES IN THE WORKING DIRECTORY, '.csv' SHOULD ONLY APPEAR AS A FILE EXTENSION !!!). This script thinks that any file with the string '.csv' in the file name is a .csv file
  dfs <- lapply(X = files, FUN = formatMicro)  # process all .csv files in the working directory into a list of lists.
  
  ## write the modified files as .csv files to an output folder
  
  # set the working directory to the folder that will collect the output files (should be empty prior to running this script)
  setwd(output.WorkingDirectory)
  
  for (j in 1:length(dfs)) {
      datavals <- dfs[[j]]$data
      in_name <- dfs[[j]]$filename
      out_name <- paste(paste(output.Prefixes[i],in_name,dfs[[j]]$date_lab,sep="_"),".csv",sep="")
      meta <- data.frame(out_filename=out_name,dfs[[j]][-1])
      write.csv(datavals, file = out_name, row.names = FALSE)
      # Writes the metadata
      if (length(list.files(".",pattern="metadata.txt")) == 0){
        write.table(meta,file="metadata.txt", sep=",",row.names = FALSE, append = FALSE) 
      }
      else{
        write.table(meta,file="metadata.txt", sep=",",row.names = FALSE, col.names = FALSE, append = TRUE)
      }
      file_n <- file_n+1
      print(paste("Completed processing file ",file_n," of ",nfiles))
  }
}

#  
