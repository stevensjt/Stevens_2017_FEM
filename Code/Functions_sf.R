##Eventually, a faster version of these functions might be possible with sf package##
##For now, documentation of sf is not sufficient to write an alternative function for "fill_holes" or "gBuffer".
##The read-in of files is much faster though:
#6 seconds; creates Multipolygon sf data frame (vs 33 seconds with read OGR)
#library(sf) #version 0.5-2
##hs_patches_sf <- 
##read_sf("../Large Files/GIS/BurnSev/Current/", layer="hs_patches")

#Finished