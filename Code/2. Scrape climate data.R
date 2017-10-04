##This code downloads climate data from the GridMet database for each fire
##Author: Jens Stevens; stevensjt@gmail.com
#Finished

#Table of contents:
#0. Load libraries
#1. Data management
#2. Get weather data from GridMet

####0. Load Libraries####
library(rgdal) #for readOGR, proj4string. version 1.2-5 (sp version 1.2-4)
library(stringr) #for str_sub. version 1.1.0
library(rgeos)# For gCentroid. version 0.3-21
library(ncdf4) #for processing functions to handle .nc files (e.g. nc_open). version 1.15
library(tidyverse) #for write_csv. version 1.1.1
source("./Code/Functions.R")

####1. Data management####
#Read in  data (without weather)
fire.list=read.csv("./Data/Derived/all_fires_ForAnalysis.csv")

#Read in spatial data (from Jens' computer):
#hs_patches=readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches") #CRS EPSG:3310, NAD83 CA Albers
#Read in spatial data from Repo (public users should do this):
hs_patches <- readRDS("./Data/hs_patches.RDS")
#Plot a fire to make sure it worked:
#plot(hs_patches[10,])

#Additional data processing

#If the end date is unknown but the start date is known, 
#Make end date 7 days after the start date 
#(bias towards period when fire was very likely to be burning).
#This line finds fires where CONT_DATE = 0; n = 24
fire.list[fire.list$CONT_DATE==0,"CONT_DATE"] <- 
  fire.list[fire.list$CONT_DATE==0,"IGNITION_DATE"] + 7
#This line finds fires where CONT_DATE has a known year but unknown month/day; n = 11
fire.list[str_sub(fire.list$CONT_DATE,-4)=="0000","CONT_DATE"]<-
  fire.list[str_sub(fire.list$CONT_DATE,-4)=="0000","IGNITION_DATE"] + 7

#If the end date is past the 30th of the month, 
#add 70 to get to the next month; n=35
fire.list[which(as.numeric(str_sub(fire.list$CONT_DATE,-2))>30),"CONT_DATE"]=
  fire.list[which(as.numeric(str_sub(fire.list$CONT_DATE,-2))>30),"CONT_DATE"] + 70 #

#If the end month is known but the end day is unknown, 
#make the end day the first day of the month; n = 7 (all from 2015 fires)
fire.list[grep("0800",fire.list$CONT_DATE), "CONT_DATE"] = fire.list[grep("0800",fire.list$CONT_DATE), "CONT_DATE"]+1
fire.list[grep("0900",fire.list$CONT_DATE), "CONT_DATE"] = fire.list[grep("0900",fire.list$CONT_DATE), "CONT_DATE"]+1

#Fires that don't have weather after this are 1984YNP, 1990 Megram, 1990 Sheep (which don't have known start dates or even months) and 2005 Crag (where start date is after end date).
#NA values for these fire in the weather column is acceptable

####2. Get weather data from GridMet ####
#Citation: Abatzoglou 2013 Int. J. Climatology
#Sample data request page: https://www.reacchpna.org/thredds/ncss/grid/MET/tmmx/tmmx_2016.nc/dataset.html
#Sample url: "https://www.reacchpna.org/thredds/ncss/MET/tmmx/tmmx_1984.nc?var=air_temperature&north=42.15&west=-124.43&east=-114.50&south=32.66&disableProjSubset=on&horizStride=1&time_start=1984-07-18T00%3A00%3A00Z&time_end=1984-07-25T00%3A00%3A00Z&timeStride=1&accept=netcdf"
#To determine full variable, go to data source, select a year of variable of interest, and click "NetcdfSubset" to study example link like what is above.
#Data source: https://www.reacchpna.org/thredds/reacch_climate_MET_catalog.html

#fire.list <- read.csv('./Data/Derived/all_fires_ForAnalysis_weather.csv') #Optional, if you've already run the full dataset and are revising or running for a subset.

fires.to.sample=as.character(fire.list$VB_ID) #N=483
#Bad fires: 251, 252, run these individually
for(f in c(1:250)){ #Fire for loop
  
  hs_fire <- #Read in the perimeter for a given high-severity patch.
  hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
  
  if(nrow(hs_fire@data)>0){ #CHECK Shapefile: 
    #If there is a shapefile in the database, proceed.
    hs_fire_ll <-  #Reproject to LL to use as NOAA argument in data call
      spTransform(hs_fire,CRS("+proj=longlat +datum=WGS84"))
    #Get dates and convert to proper format using regex
    sd=fire.list[f,"IGNITION_DATE"]
    sd=gsub('^(.{4})(.*)$', '\\1-\\2', sd) #Insert dash (crazy regex)
    sd=gsub('^(.{7})(.*)$', '\\1-\\2', sd) #Insert dash
    ed=paste(fire.list[f,"CONT_DATE"])
    ed=gsub('^(.{4})(.*)$', '\\1-\\2', ed) #Insert dash (crazy regex)
    ed=gsub('^(.{7})(.*)$', '\\1-\\2', ed) #Insert dash
    
    vs <- #Variables selected; High temperature, low temp, high humidity, burn index
      c("tmmx","tmmn","rmax","bi") 
    vs_full <- #API call needs full variable name as well
      c("air_temperature","air_temperature","relative_humidity","burning_index_g")
    for(index in 1:length(vs)){ #Variable for loop
      v=vs[index]
      v_full=vs_full[index]
      
      #Key function: Create a string with the url for data download
      #This string includes the specific area (lat/long) and time window of interest. 
      #The area is specified by a bounding box,
      #with a 0.2 degree buffer around the fire centroid.
      #The resolution of the data is 0.04 degrees per pixel-width.
      v_link=paste0("https://www.reacchpna.org/thredds/ncss/MET/",v,
                    "/",v,"_",fire.list[f,"FIRE_YEAR"],".nc?var=",v_full,
                    "&north=",gCentroid(hs_fire_ll)$y+0.2,
                    "&west=",gCentroid(hs_fire_ll)$x-0.2,
                    "&east=",gCentroid(hs_fire_ll)$x+0.2,
                    "&south=",gCentroid(hs_fire_ll)$y-0.2,
                    "&disableProjSubset=on&horizStride=1",
                    "&time_start=",sd,"T00%3A00%3A00Z",
                    "&time_end=",ed,"T00%3A00%3A00Z&timeStride=1&accept=netcdf")
      
      dest <-  paste0("./Data/climate_nc/",v,".nc" )
      #Optional: Add lat/long centroid to data frame to check accuracy.
      #fire.list[f,"lat"]=gCentroid(hs_fire_ll)$y
      #fire.list[f,"long"]=gCentroid(hs_fire_ll)$x
      tmp <- #Download the relevant nc file from GridMet for the appropriate point and time
        try(download.file(url=v_link,destfile=dest), silent=T)
      if(class(tmp)!="try-error"){ #CHECK download error: 
        #if the download produced an error (e.g. fires 251, 252), 
        #you will have to do it manually so skip the next bit.
        #otherwise, extract climate maxima/minima from the downloaded file.
        ncin <- nc_open(dest) #Open the file you just downloaded
        lat <- ncvar_get(ncin,"lat",verbose=F) #Several pixels in sample
        lat_target <- #Find closest latitude pixel to fire centroid
          which.min(abs(lat - gCentroid(hs_fire_ll)$y)) 
        lon <- ncvar_get(ncin,"lon",verbose=F)
        lon_target <- #Find closest latitude to fire centroid
          which.min(abs(lon - gCentroid(hs_fire_ll)$x)) 
        v_array <- ncvar_get(ncin,v_full)
        if(class(v_array)=="matrix"){ 
          #If there was only one burn day, so there's a matrix instead of an array
          #Duplicate the arrayso there's no error produced in maximum calculation
          v_array=replicate(2,v_array,simplify="array")
        }
        if(v=="tmmx"){
          #Get maximum high temperature during the burn window
          fire.list[f,paste0("max_",v)] <- #Convert from K to C
            max(v_array[lat_target,lon_target,])-273.15 
        }
        if(v=="tmmn"){
          #Get maximum low temperature during the burn window
          fire.list[f,paste0("max_",v)] <- #Convert from K to C
            max(v_array[lat_target,lon_target,])-273.15
        }
        #Get minimum high RH during the burn window
        if(v=="rmax"){
          fire.list[f,paste0("min_",v)] <-
            min(v_array[lat_target,lon_target,]) 
        }
        if(v=="bi"){
          #Get max burn index during burn window
          fire.list[f,paste0("max_",v)] <-
            max(v_array[lat_target,lon_target,]) 
        }
      } #END CHECK download error
    } #END variable for loop
    gc()
    
  } #END CHECK shapefile
  print(f)
  gc()
} #END fire for loop

#Deer fire is wierd; very small perimeter and has negative BI, so set to NA
fire.list[which(fire.list$max_bi<0),"max_bi"]=NA 

#write_csv(fire.list,'./Data/Derived/all_fires_ForAnalysis_weather.csv') #Save file