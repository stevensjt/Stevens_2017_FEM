##This code reads in mapped stand-replacing patches for a specified set of fires
##This code then calculates SDC for those patches
##Author: Jens Stevens; stevensjt@gmail.com
#Finished

#Table of contents:
#0. Load libraries
#1. Load high severity patches and list of fires to analyze
#2. Run Internal Buffering ['decay()'], create "Long" dataset
#3. Calculate SDC 
#4. Save datasets

####0. load libraries####
library(rgdal) #for readOGR. version 1.2-5 (sp version 1.2-4)
library(tidyverse) #for ggplot + tidyverse functions. version 1.1.1
require(rgeos) #for createSPComment. version 0.3-21

####1. Load high severity patches and list of fires to analyze####

##1a. Load high-severity patches
#The input "hs_patches.shp" that is read in here contains all fires from the USFS database, available at: https://www.fs.usda.gov/detail/r5/landmanagement/gis/?cid=stelprd3804878

#The process of extracting the high severity patches from the raw USFS severity shapefiles in R is very time consuming. I'm doing this in ArcGIS (on PC) and creating shapefile called "hs_patches". The high-severity patches are created by setting the RdNBR value associated with 90% basal area mortality as the minimum threshold value.

#Reading data from Jens' computer:
#hs_patches <- #33 seconds; creates Large SpatialPolygonsDataFrame
#readOGR("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")
#gc() #Clear space
#Write the shapefile as a compressed .RDS file, save to repo directory.
#saveRDS(hs_patches,"./Data/hs_patches.RDS")

#Reading data from Repo (public users should do this):
hs_patches <- readRDS("./Data/hs_patches.RDS")
#Plot a fire to make sure it worked:
#plot(hs_patches[10,])

#1b. Read in raw fire list file from Jay Miller
fire.list <- 
  read.csv("./Data/Raw/fires_usfs.csv") %>%
  filter(Veg == "Forest",(PCT_FS > 0.5 | (PCT_FS < 0.5 & AGENCY == "NPS") ), 
                          FIRE_YEAR < 2016 ) 
fire.list <- fire.list[fire.list$VB_ID %in% hs_patches$VB_ID,]
#This dataset now has fires >80 ha from 1984 through 2015 that burned predominantly through conifer forest vegetation, and were either at least 50% on Forest Service lands or were at least 50% on Park Service lands (managed by NPS), and have a corresponding shapefile of high-severity patches. Missing 10 fires from 2016 that haven't been mapped yet; because 2016 was an incomplete year we excluded it from analyses.
#N = 477 fires


##1c. Create character vector of fires to sample
fires.to.sample <- as.character(fire.list$VB_ID)[1:10] #Sample all fires
#fires.to.sample <- c("1987EAST","2008CARIBOU") #Sample specific fires

####2. Run Internal Buffering ['decay()'], create "Long" dataset####
Sys.time() #Takes ~11 hours in parallel (Rim takes 2 hours). 
#Speed up by changing buf_inc to 20 or 30; SDC is about the same.

for(f in c(1:length(fires.to.sample))){
  hs_fire <- #Grab the fth fire in the list, load in the spatial layer 
    hs_patches[hs_patches$VB_ID==fires.to.sample[f],]
  hs_fire <- # Fill holes <0.81 ha (from Functions.R)
    fill_holes(hs_fire=hs_fire)
  hs_fire2 <-   #This is needed to buffer from internal "holes" within the patch > 0.81 ha.
    createSPComment(hs_fire)
  #plot(hs_fire,col="darkred",border="transparent")
  Sys.time() #CHECKME 33 sec for East Fire.
  decay.table <- # internal buffering at 10 m increments (from Functions.R)
    decay(hs_fire=hs_fire2,buf_max=1000,buf_inc=10,name=fires.to.sample[f])
  Sys.time()
  ifelse(f==1, fires_long <- decay.table, fires_long <- rbind(fires_long,decay.table) )
  print(paste(fires.to.sample[f],Sys.time()))
  gc()
}

####3. Calculate SDC####
#Calculate SDC for each fire (fast), and summarize to a single row per fire 
fires_long <- calculate.sdc(fires_long)
summary_fires <- 
  fires_long %>%
  group_by(name) %>%
  summarise(sdc = mean(sdc))
names(summary_fires)[1] <- "VB_ID"
#Add the sdc parameter to the fire.list file, with a single row per fire 
#This will get saved as the "ForAnalysis" dataset.
fire.list <- merge(fire.list,summary_fires[,c("VB_ID","sdc")])

####4.Save datasets####
#Save as fires_long as "Long" dataset, including the sdc values 
#(use this file for plotting curves).
write_csv(fires_long,"./Data/Derived/all_fires_Long.csv")

#Save fires.list as "For Analysis" short dataset, one row per fire with SDC.
write_csv(fire.list,path="./Data/Derived/all_fires_ForAnalysis.csv")