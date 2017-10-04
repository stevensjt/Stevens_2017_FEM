##This code creates and plots circles with different SDC values for Appendix Figure A2
##Author: Jens Stevens; stevensjt@gmail.com
#Finished

#Table of contents:
#0. Load Libraries
#1. Create and plot circles with different SDC values

####0. Load Libraries####
library(gridExtra) #for grid.arrange. version 2.2.1
library(rgeos) # For gBuffer. version 0.3-21

####1. Create and plot circles with different SDC values####
c.df <- data.frame(c.rad=c(20,40,60,80,100,150,200,250,300,400,500,1000),
                   sdc=NA)
c.polys=c.rad=c.plots=list()
for(r in 1:nrow(c.df)){
  c.rad[[r]] <- c.df[r,"c.rad"]
  coords_sp <- SpatialPoints(data.frame(x=0,y=0),CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")) #CRS EPSG:3310, NAD83 CA Albers
  c.poly <- gBuffer(coords_sp,width=c.rad[[r]],quadsegs = 20)
  c.poly.decay <- decay(c.poly,buf_inc=10,buf_max=1000,name=as.character(c.rad[[r]]))
  c.df[r,"sdc"] <- calculate.sdc(c.poly.decay)[,"sdc.name"] %>% unique()
  c.polys[[r]] <- c.poly
  c.plots[[r]] <-
    ggplot()+
    geom_path(data=c.poly,aes(x=long,y=lat))+
    geom_segment(aes_string(x=0,y=0,xend=c.df[r,"c.rad"],yend=0))+
    annotate("text",x=-1000,y=900,label=
               paste("r =",c.rad[[r]]),hjust=0,vjust=0,size=3)+
    annotate("text",x=-1000,y=-1000,label=
               paste("a =\n",round(0.0001*pi*c.rad[[r]]^2,2), "ha"),hjust=0,vjust=0,size=3)+
    lims(x=c(-1000,1000),y=c(-1000,1000))+
    labs(title=paste0("sdc = ",c.df[r,"sdc"], 
                      "\n ln(sdc) = ", round(log(as.double(c.df[r,"sdc"])),3)))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
}

png(file = paste0("./Figures/FigA2_",Sys.Date(),".png"),width=8,height=12,units="in",res=200)
do.call("grid.arrange", c(c.plots, ncol=3))
dev.off()