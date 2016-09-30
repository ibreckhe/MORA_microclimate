####Script to download MesoWest (RAWS,SNOTEL,NWS) station data for the Mt. Rainier vicinity
##using the MesoWest API.

library(jsonlite)
library(RCurl)


stationMeta <- getURLContent("http://api.mesowest.net/v2/stations/metadata?&token=2c7ec72169f2461ba048989f94e1e11d&complete=1&radius=46.853634,%20-121.764424,%2024.8548")
stationMeta_df <- fromJSON(stationMeta)$STATION[,-3]
write.csv(stationMeta_df,"MW_Station_metadata_40km.csv",row.names=FALSE)
stationIDs <- stationMeta_df$STID

##Function for grabbing the data.
get_mesowest_data <- function(api_token,date_start,date_end,var,stationID,filepath,filename){
  api_list <- list(start="http://api.mesowest.net/v2/stations/timeseries?&token=",
                   token=paste(api_token,"&start=",sep=""),
                   start=paste(date_start,"&end=",sep=""),
                   end=paste(date_end,"",sep=""),
                   units="&units=metric",
                   stid=paste("&stid=",stationID,sep=""),
                   var=paste("&vars=",var,sep=""))
  api_string <- paste(api_list,collapse="")
  
  ##calls the API:
  result <-getURLContent(api_string)
  result_str <- gsub("\\{.*\\}","",result)
  
  ##Writes file to disk if the call returns data and checks the formatting of the resulting file.
  fullpath <- paste(filepath,filename,sep="")
  if(nchar(result_str)==0){
    print(paste("No data returned for station",stationID))
  }else{
    writeLines(result_str,fullpath)
    ts_data <- data.table::fread(fullpath,sep=",")
    if(nrow(ts_data)>3 & colnames(ts_data)[1]=="Station_ID"){
      print(paste("Successfully read",nrow(ts_data),"observations for station",stationID))
    }else{
      print(paste("Few (<4) observations for station ",stationID))
    }
  }
}

##Loops through all the stations and grabs temp data.
api_token <- "5e9f0a9379554e1d950b99f4c2f5a087"
date_start <- 200901011201
date_end <- 201512311201

##Gets air temp.
for (i in 1:length(stationIDs)){
  flush.console()
  paste("Downloading data for station",i,"of",length(stationIDs))
  var<-"air_temperature"
  filepath <- "~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/MesoWest/"
  filename <- paste(stationIDs[i],"_",var,"_",date_start,"_",date_end,".csv",sep="")
  get_mesowest_data(api_token=api_token,date_start=date_start,
                    date_end=date_end,var="air_temp",stationID=stationIDs[i],
                    filepath=filepath,filename=filename)
}
##Gets relative humidity
for (i in 1:length(stationIDs)){
  flush.console()
  paste("Downloading data for station",i,"of",length(stationIDs))
  var<-"relative_humidity"
  filepath <- "~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/MesoWest/"
  filename <- paste(stationIDs[i],"_",var,"_",date_start,"_",date_end,".csv",sep="")
  get_mesowest_data(api_token=api_token,date_start=date_start,
                    date_end=date_end,var=var,stationID=stationIDs[i],
                    filepath=filepath,filename=filename)
}

##Gets wind speed.
for (i in 1:length(stationIDs)){
  flush.console()
  paste("Downloading data for station",i,"of",length(stationIDs))
  var<-"wind_speed"
  filepath <- "~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/MesoWest/"
  filename <- paste(stationIDs[i],"_",var,"_",date_start,"_",date_end,".csv",sep="")
  get_mesowest_data(api_token=api_token,date_start=date_start,
                    date_end=date_end,var=var,stationID=stationIDs[i],
                    filepath=filepath,filename=filename)
}

##Gets Solar Radiation.
for (i in 1:length(stationIDs)){
  flush.console()
  paste("Downloading data for station",i,"of",length(stationIDs))
  var<-"solar_radiation"
  filepath <- "~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/MesoWest/"
  filename <- paste(stationIDs[i],"_",var,"_",date_start,"_",date_end,".csv",sep="")
  get_mesowest_data(api_token=api_token,date_start=date_start,
                    date_end=date_end,var=var,stationID=stationIDs[i],
                    filepath=filepath,filename=filename)
}

##Gets Snow Water Equivalent.
for (i in 1:length(stationIDs)){
  flush.console()
  paste("Downloading data for station",i,"of",length(stationIDs))
  var<-"snow_water_equiv"
  filepath <- "~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/MesoWest/"
  filename <- paste(stationIDs[i],"_",var,"_",date_start,"_",date_end,".csv",sep="")
  get_mesowest_data(api_token=api_token,date_start=date_start,
                    date_end=date_end,var=var,stationID=stationIDs[i],
                    filepath=filepath,filename=filename)
}

##Gets pressure.
for (i in 1:length(stationIDs)){
  flush.console()
  paste("Downloading data for station",i,"of",length(stationIDs))
  var<-"pressure"
  filepath <- "~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Microclimate/raw/MesoWest/"
  filename <- paste(stationIDs[i],"_",var,"_",date_start,"_",date_end,".csv",sep="")
  get_mesowest_data(api_token=api_token,date_start=date_start,
                    date_end=date_end,var=var,stationID=stationIDs[i],
                    filepath=filepath,filename=filename)
}



                 