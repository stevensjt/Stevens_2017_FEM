##This code analyzes trends in SDC, 
##looks at variables that explain variation in SDC, 
##and examines individual example fires
##Author: Jens Stevens; stevensjt@gmail.com
#Finished

#Table of contents:
#0. Load Libraries
#1. Read in, process, and filter data for analysis
#2. Exploratory analyses
#3a. SDC model selection
#3b. Regression tree
#4a. Trends over time
#4b. Trends over time by region
#5a. Identify specific example fires
#5b. Plot specific example fires
#6. Plot SDC vs class, agency, percent high-severity, and fire size
#7. Analyze "forest loss" by agency

####0. Load Libraries####
library(tidyverse) #for %>%, write_csv, etc. version 1.1.1
library(ReporteRs) #for cellProperties etc. version 0.8.8
library(glmulti) #for glmulti. version 1.0.7
library(rpart) #for rpart. version 4.1-10
library(rpart.plot) #for prp. version 2.1.2
library(RColorBrewer) #For brewer.pal. version 1.1-2
library(grid) #for viewport. version 3.3.2
library(gridExtra) #for grid.arrange. version 2.2.1
library(car) #for durbinWatsonTest. version 2.1-4
library(rgdal) #for readOGR. version 1.2-5 (sp version 1.2-4)

source("./Code/Functions.R") #Need for fill_holes, to create Fig. 6.

####1. Read in, process, and filter data for analysis####
#This section can be skipped; the core data product sdc_data.csv produced by this step can be read in during step 2.
#d <- read.csv("./Data/Derived/all_fires_ForAnalysis_weather.csv")

d <- d[-which(d$FIRE_YEAR==2016),] #N=477

#modify variable names for plotting purposes
names(d) <- c("VB_ID","ID_Num","fire_name","fire_year","class","veg_type",
              "agency","ICS_code","pct_fs", #<-MGMT variables
              "firesize_ha","BA90_ha","BA90_pct","max_patch_ha",
              "ignition_date","contain_date", #<-FS stats
              "sdc","max_tmmx","max_tmmn",
              "min_rmax","max_bi") #<-spatial stats & weather

#Rename agency and class codes for clarity
d$agency=as.character(d$agency)
d[d$agency%in%c("BIA","CCO"), "agency"]= NA #Get rid of a few underrepresented agencies
d[grep("USF",d$agency), "agency"]= "USFS" #Recode "USF" as "USFS" for consistency in text; this also classifies five fires with agency of "USF/NPS" as "USFS" because they had a National Forest as the ICS code (Jay's advice).
d$agency <- factor(d$agency)
d$class <- factor(d$class,levels=c("no","yes"), labels =c("SUP","WFU"))

#Add region variable (Sierra Nevada/Northwest) based on ICS code.
d$region <- ifelse(d$ICS_code%in%c("KNF","MEU","MNF","SHF","SHU","SRF"),"NW","SCSN")
d$region <- factor(d$region,levels=c("SCSN","NW"))

#Two funky fires that had 100 percent humidity and very low temps; suspect potential bias in extraction of meteorological data so excluding met data for these fires.
d[which(d$min_rmax==100),c(17:20)]=NA 

#If running for the first time, save core data frame to include "region" variable.
#write_csv(d,"./Data/sdc_data.csv") 

####2. Exploratory analyses####
d <- read.csv("./Data/sdc_data.csv")
#2a: Check for normality
#hist(d$sdc)
#hist(log(d$sdc)) #note in R: log = ln

#2b: Summary table (Table 2)
Table2 <- 
  d %>%
  group_by(agency,class) %>%
  summarize(N = length(sdc), min_ha = min(firesize_ha), median_ha = median(firesize_ha), 
            max_ha = max(firesize_ha,na.rm=T), median_fireyear = round(median(fire_year,na.rm=T),0), 
            mean_max_tmmx = mean(max_tmmx,na.rm=T), mean_max_bi = mean(max_bi,na.rm=T), 
            mean_max_tmmn = mean(max_tmmn,na.rm=T), mean_min_rmax = mean(min_rmax,na.rm=T)
  )

Table2[,c(8:11)] <- round(Table2[,c(8:11)],1)
Table2[,"agency"] <- as.character(Table2[,"agency"][[1]])
names(Table2) <- c("agency","class","N","min size \n(ha)","median \nsize (ha)","max size \n(ha)","median \nfire year",
                   "mean maximum \nhigh temperature","mean maximum\nburn index","mean maximum \nlow temperature", 
                   "mean minimum \nhigh humidity")
Table2[is.na(Table2$agency),"agency"] <- "NA"
# Set up general table properties and formatting
cell_p = cellProperties(padding.right=3, padding.left=3)
par_p = parProperties(text.align="right")
# Make Table
ft = FlexTable(Table2, header.columns=FALSE, body.cell.props=cell_p, body.par.props=par_p)
ft = addHeaderRow(ft, text.properties=textBold(), names(Table2),
                   par.properties=parCenter())
ft #Save as HTML Table 2

####3a: SDC model selection####
#tmax and tmin are correlated, so just using tmax
#WFU and AGENCY are clearly the most important variables
#Best AIC: FIRE_YEAR, AGENCY, WFU, max_tmmx, max_tmmn (890.12). tmmx and tmmn are correlated, so they offset each other in the model but there's a slightly greater effect size for tmmx
#Equivocal model removes tmmn (and reduce the effect size of tmmx (890.73));
#tmax and tmin are correlated, so just using tmax is acceptable
#Equivocal model replaces FIRE_YEAR with burn index (they are correlated) (890.86)
#Big jumps happen from 26:27 (can't really explain; adding FIRE_YEAR to an equation that had all the weather variables except tmax) and 30:31 (can't really explain; adding tmmn and rmax to equation that just had FIRE_YEAR plus WFU and AGENCY). Generally it's a pretty gradual shift.

potential_parms=c("fire_year","class","agency","region","max_tmmx","max_tmmn","min_rmax","max_bi")
m2c <- glmulti(y="log(sdc)", 
        xr=potential_parms,
        data=d,level=1,method="h") #If running interactions (level 2), try genetic algorithm (method = "g")

summary(m2c@objects[[2]])

#Make table 1: Model summary
m_max <- 5
x=matrix(NA,nrow=length(c("AIC",rev(rownames(coef(m2c) ) ) ) ),ncol=m_max+1)
x[,1] <-c("AIC",sort(rownames(coef(m2c))))

for(m in c(1:m_max)){
  c <- round( summary(m2c@objects[[m]])$coefficients[,1]
              [order(names(summary(m2c@objects[[m]])$coefficients[,1]))], 3)
  x[which(x[,1]%in%names(c)),m+1] <- c
  x[1,m+1] <- round(AIC(m2c@objects[[m]]),3)
}
x = x[c(1,2,3,4,5,6,9,8,11,7,10),] #Order weather variables by their inclusion in the model
# Set up general table properties and formatting
cell_p = cellProperties(padding.right=3, padding.left=3)
par_p = parProperties(text.align="right")
# Create table
ft = FlexTable(x, header.columns=FALSE, body.cell.props=cell_p, body.par.props=par_p)
ft = addHeaderRow(ft, text.properties=textBold(), c("","Model #"),
                  colspan=c(1,m_max), par.properties=parCenter())
ft = addHeaderRow(ft, text.properties=textBold(), c("Model AIC \n /coefficients",c(1:m_max)),
                  colspan=rep(1,times=m_max+1), par.properties=parCenter())
ft #Save as HTML Table 1


####3b: Regression tree with best model from above####
#Choosing the simplest (2nd best) model; within one AIC point of the best model. 
#http://blog.revolutionanalytics.com/2013/06/plotting-classification-and-regression-trees-with-plotrpart.html
#Recode variables for regression tree viz
rtree.d <- d
names(rtree.d)[which(names(rtree.d)=="max_tmmx")] <- "max_high_temp"
levels(rtree.d$class) <- c("suppression","wildland fire use")
#formula(m2c@objects[[2]]) #Get formula and modify max_high_temp as below
tree.1 <- rpart(formula = "log(sdc) ~ 1 + class + agency + fire_year + max_high_temp", data=rtree.d)

split.labs <- function(x, labs, digits, varlen, faclen) {
  sapply(labs, function(lab) 
    if (grepl("fire_year", lab)) {
      rhs <- sub(".* ", "", lab);
      lab <- sub(rhs, ceiling(as.numeric(rhs))+1, lab) #+1 here is the modification
    } else lab)
} #Note: The default printing takes the fire_year split at 2010.5 (see print(tree.1)) and incorrectly rounds it to >=2010 instead of >=2011. The "+1" modification above is the only way I could figure out to make this round correctly.

png(file = paste0("./Figures/Fig2_",Sys.Date(),".png"),width=4.5,height=4,units="in",res=500)

h <- #Histogram of all fires to put ln(SDC) in context
  ggplot(d)+
  geom_histogram(aes(log(sdc)), binwidth = 0.333,
                 fill=c("darkred","darkred",brewer.pal(10,"RdYlBu")),
                 col="black")+
  geom_vline(xintercept=c(-3.8,-5.1),col="black", lty=2)+
  labs(x = "ln(sdc)", title = "all fires")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7))

prp(tree.1, varlen = 0, faclen = 0, type = 3, extra = 1,
    box.col = brewer.pal(10,"RdYlBu")[c(4,4,4,4,4,4,4,5,5,7,6,6,6,8,7)],
    #box.col is custom to match the histogram; values match tree.1$frame$yval
    #extra = 1, under = TRUE,
    #yesno.yshift=-1, #Doesn't apply if using type = 3
    clip.right.labs = FALSE, #Only applies if using type = 3
    mar = c(2,2,1,2), 
    split.fun = split.labs, main="                            ln(SDC)") # Plot the tree
#varlen: no limit to variable name length
#faclen: no limit to factor name length
#Type: How to place the labels (good types = 0,3)
#yesno.yshift: move yes/no labels further down

print(h, vp=viewport(.8, .25, .4, .45))
dev.off()


#Investigate the wierd max high temp > 39 result. It can be explained by the fact that all (18) of these fires were in the northwest, where complex topography can give more complex stand-replacing fire dynamics, and most (10) were in 1987 which was a warm year but where topography might have regulated stand-replacing effects.
#d_tmp <- d[d$class=="SUP" & d$max_tmmx >= 24 & d$fire_year<2011 & d$max_tmmx >= 39 & !is.na(d$fire_name),]
#d_tmp <- d_tmp[!is.na(d_tmp$fire_name),] 

####4a: Trends over time####

d.annual <-
  group_by(d,fire_year) %>%
  summarise(sdc=mean(sdc), log_sdc=mean(log(sdc)), BA90_pct = mean(BA90_pct), 
            max_bi=mean(max_bi,na.rm=T),max_tmmx=mean(max_tmmx,na.rm=T))

sdc <- 
  ggplot(d.annual,aes(x=fire_year,y=log_sdc))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x=" ", y= "ln(sdc)",title="a")+
  annotate("text", x=2015, y=-3.4, size = 4, label = paste("R^2 == ", 0.11), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=-3.7, size = 4, label = paste("P = ", 0.058),hjust=1)+
  theme_bw()+
  theme(plot.title=element_text(size=14))
bi <- 
  ggplot(d.annual,aes(x=fire_year,y=max_bi))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x=" ", y= "burn index", title = "b")+
  annotate("text", x=2015, y=55, size = 4, label = paste("R^2 == ", 0.32), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=52, size = 4, label = paste("P = ", 0.001),hjust=1)+
  theme_bw()
tmmx <- 
  ggplot(d.annual,aes(x=fire_year,y=max_tmmx))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "maximum temperature (C)", title = "c")+
  annotate("text", x=2015, y=26, size = 4, label = paste("R^2 == ", 0.10), parse = TRUE,hjust=1)+
  annotate("text", x=2015, y=24, size = 4, label = paste("P = ", 0.077),hjust=1)+
  theme_bw()

png(file = paste0("./Figures/Fig4_",Sys.Date(),".png"),width=3.5,height=9,units="in",res=500)
grid.arrange(sdc,bi,tmmx,ncol=1)
dev.off()

#Changes over time in log(sdc) marginally significant
m <- lm(log(sdc)~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m) #check for temporal autocorrelation (none found)

#Changes over time in burn index - significantly more positive
m <- lm(max_bi~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m) #check for temporal autocorrelation (none found)

#Changes over time in max temperature marginally significant
m <- lm(max_tmmx~fire_year,data=d.annual)
summary(m)
durbinWatsonTest(m) #check for temporal autocorrelation (none found)

####4b: Trends over time, subset by region####
d.annual.region <-
  group_by(d,fire_year, region) %>% 
  summarise(sdc=mean(sdc), log_sdc=mean(log(sdc)),  BA90_pct = mean(BA90_pct), 
            max_bi=mean(max_bi,na.rm=T),max_tmmx=mean(max_tmmx,na.rm=T)) 

d.annual.region$region <- factor(d.annual.region$region,levels=c("SCSN","NW"))

sdc_region <- 
  ggplot(d.annual.region,aes(x=fire_year,y=log_sdc,col=region))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="fire year", y= "ln(sdc)")+
#  annotate("text", x=2015, y=-3.4, size = 4, label = paste("R^2 == ", 0.14), parse = TRUE,hjust=1)+
#  annotate("text", x=2015, y=-3.7, size = 4, label = paste("P = ", 0.047),hjust=1)+
  theme_bw()+
  theme(plot.title=element_text(size=14),legend.position = "none")

png(file = paste0("./Figures/FigA3_",Sys.Date(),".png"),width=3,height=3,units="in",res=500)
sdc_region
dev.off()

summary(lm(log_sdc~fire_year,data=d.annual.region[d.annual.region$region=="NW",]))
summary(lm(log_sdc~fire_year,data=d.annual.region[d.annual.region$region=="SCSN",]))
#significant trend for Sierra Nevada but not northwestern CA.

####5a. Identify specific example fires####
#OK to skip 5a; example fires identified here and referenced by name in 5b.
#Run the chunk below for BA_range <- c(0,10), c(10,20), c(20,30), c(30-40), and c(45-55). 
BA_range <- c(0,10) #Need to change this manually to each of the BA_range values above.
BA_name <- paste("BA",BA_range[1],BA_range[2],"candidate",sep="_")
for(r in 1:nrow(d)){ #Find fires with BA90 in a desired range
  if(between(d$BA90_pct[r],BA_range[1],BA_range[2])){
    d[r,BA_name] = TRUE
  } else d[r,BA_name] = FALSE
}
FS_range <- quantile(as.vector(na.exclude(d[d[,BA_name] == TRUE,"firesize_ha"])),c(0.45,0.55))
for(r in 1:nrow(d)){ #Identify those fires of similar size
  if(d[r,BA_name] & between(d$firesize_ha[r],FS_range[1],FS_range[2])){
    d[r,BA_name] = TRUE
  } else {d[r,BA_name] = FALSE}
}
for(r in 1:nrow(d)){ #Identify those fires of similar size
  if(d[r,BA_name] & between(d$firesize_ha[r],FS_range[1],FS_range[2])){
    d[r,BA_name] = TRUE
  } else {d[r,BA_name] = FALSE}
}
d[d$sdc %in% c( max(d[d[,BA_name]==TRUE,"sdc"]), min(d[d[,BA_name]==TRUE,"sdc"]) ),"Plot_candidates"]=TRUE

####5b.Plot  specific example fires####
#Specify the example fires from those identified in 5a.
d[as.character(d$VB_ID) %in% 
    c("2001HIGHWAY","1990RECER","1999HORTON2","2009GOOSE","1989RACK","2006BOULDER_CMPLX",
      "2001STREAM","2004SIMS"),"Plot_candidates"] <-  c(5,2,3,1,7,8,6,4)
d[which(is.na(d$Plot_candidates)),"Plot_candidates"] <- ""


hs_patches <- readRDS("./Data/hs_patches.RDS") #CRS EPSG:3310, NAD83 CA Albers; Repo users can access this.
full_shapes <- readOGR("../Large Files/GIS/BurnSev/Current/", layer="SampleFiresFullPerims") #CRS EPSG:3310, NAD83 CA Albers
#Using full_shapes above on reviewer recommendation to display high severity areas within context of full fire perimeter. Dataset not available in public repo.

fires.to.plot.names=as.character(d[ order(d[,"Plot_candidates"])[470:477] , "VB_ID"])
pd.a <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[1],])
p.a=
  ggplot()+
  geom_polygon(data=pd.a,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[1],]),
               aes(x=long-min(pd.a$long),y=lat-min(pd.a$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="1. Highway Fire (2001)",x=" ",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[1],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[1],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.b <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[2],])
p.b=ggplot()+
  geom_polygon(data=pd.b,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[2],]),
               aes(x=long-min(pd.b$long),y=lat-min(pd.b$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="2. Recer Fire (1990)",x=" ",y=" ") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[2],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[2],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.c <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[3],])
p.c=ggplot()+
  geom_polygon(data=pd.c,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[3],]),
               aes(x=long-min(pd.c$long),y=lat-min(pd.c$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="3. Horton Fire (1999)",x=" ",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[3],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[3],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.d <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[4],])
p.d=ggplot()+
  geom_polygon(data=pd.d,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[4],]),
               aes(x=long-min(pd.d$long),y=lat-min(pd.d$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="4. Goose Fire (2009)",x=" ",y=" ") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[4],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[4],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.e <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[5],])
p.e=ggplot()+
  geom_polygon(data=pd.e,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[5],]),
               aes(x=long-min(pd.e$long),y=lat-min(pd.e$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="5. Rack Fire (1989)",x=" ",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[5],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[5],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.f <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[6],])
p.f=ggplot()+
  geom_polygon(data=pd.f,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[6],]),
               aes(x=long-min(pd.f$long),y=lat-min(pd.f$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="6. Boulder Complex (2006)",x=" ",y=" ") +
  annotate("text",x=0,y=4000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[6],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[6],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.g <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[7],])
p.g=ggplot()+
  geom_polygon(data=pd.g,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[7],]),
               aes(x=long-min(pd.g$long),y=lat-min(pd.g$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="7. Stream Fire (2001)",x="meters",y="meters") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[7],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[7],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

pd.h <- fortify(full_shapes[full_shapes$VB_ID==fires.to.plot.names[8],])
p.h=ggplot()+
  geom_polygon(data=pd.h,
               aes(x=long-min(long),y=lat-min(lat),group=group),fill="gray",col="gray")+
  geom_polygon(data=fill_holes(hs_patches[hs_patches$VB_ID==fires.to.plot.names[8],]),
               aes(x=long-min(pd.h$long),y=lat-min(pd.h$lat),group=group),col='darkred',fill='darkred')+
  xlim(0,10000) + ylim(0,10000) + coord_fixed()+
  labs(title="8. Sims Fire (2004)",x="meters",y=" ") +
  annotate("text",x=0,y=8000,label=
             paste("sdc =",round(d[d$VB_ID==fires.to.plot.names[8],"sdc"],4),
                   "\nln(sdc) = ",round(log(d[d$VB_ID==fires.to.plot.names[8],"sdc"]),2) ),
           hjust=0,vjust=0,size=3)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12,angle = 0),
        axis.title=element_text(size=13),
        plot.title=element_text(size=14, hjust=0.5))

png(file = paste0("./Figures/Fig6_",Sys.Date(),".png"),width=6,height=10,units="in",res=500)
grid.arrange(p.a,p.b,p.c,p.d,p.e,p.f,p.g,p.h,ncol=2)
dev.off()

####6: Plot SDC vs class, agency, pct_hs, and fire size#### 
#Make sure you've run #5b first.

#Trouble with plotting NA's as gray
#df <- data.frame (V1=factor(c("A","B","A","B",NA)),x=c(1:5),y=c(1:5))
#ggplot(df,aes(x=x,y=y,col=V1))+
#  geom_point()

class_pct <-
  ggplot(na.omit(d[,c("BA90_pct","sdc","class","Plot_candidates")]),aes(x=BA90_pct,y=log(sdc),col=class))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_text(aes(label=Plot_candidates),hjust=1, vjust=1,size=6,col="black",fontface="bold")+
  labs(y= "ln(sdc)",x=" ", title = "            percent high-severity\na")+
  annotate("text", x=60, y=-3, label = paste("R^2 == ", 0.67), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.position = 'none')
class_ha <-
  ggplot(na.omit(d[,c("firesize_ha","sdc","class")]),aes(x=log(firesize_ha),y=log(sdc),col=class))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y= " ", x=" ", title = "                       fire area\nb")+
  annotate("text", x=10, y=-3, label = paste("R^2 == ", 0.22), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13))
agency_pct <-
  ggplot(na.omit(d[,c("BA90_pct","sdc","agency")]),aes(x=BA90_pct,y=log(sdc),col=agency))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"),guide=FALSE)+
  geom_smooth(method="lm")+
  labs(y= "ln(sdc)", x="% high-severity", title = "c")+
  annotate("text", x=60, y=-3, label = paste("R^2 == ", 0.67), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13))
agency_ha <-
  ggplot(na.omit(d[,c("firesize_ha","sdc","agency")]),aes(x=log(firesize_ha),y=log(sdc),col=agency))+
  geom_point()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  geom_smooth(method="lm")+
  labs(y= " ", x= "ln(area [ha])", title = "d")+
  annotate("text", x=10, y=-3, label = paste("R^2 == ", 0.21), parse = TRUE) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13))

#Test joint effects of (% high severity or fire size) and (agency or class) on sdc variation
summary(lm(log(sdc)~BA90_pct+class,data=d))
summary(lm(log(sdc)~firesize_ha+class,data=d))
summary(lm(log(sdc)~BA90_pct+agency,data=d))
summary(lm(log(sdc)~BA90_pct+agency,data=within(d, agency <- relevel(agency, ref = 3) ) ) ) #change factor order
summary(lm(log(sdc)~firesize_ha+agency,data=d))
summary(lm(log(sdc)~firesize_ha+agency,data=within(d, agency <- relevel(agency, ref = 3) ) ) ) #change factor order

png(file = paste0("./Figures/Fig3_",Sys.Date(),".png"),width=7.5,height=10,units="in",res=500)
grid.arrange(class_pct,class_ha,agency_pct,agency_ha,ncol=2,widths=c(0.44,0.56))
dev.off()

####7: Forest loss by agency####
#Calculate core patch area (120 m buffer in from edge)
d$P <- 1/(10^(d$sdc*120))
d$CPA_120 <- d$BA90_ha * d$P

d_CPA <-
  d %>%
  group_by(fire_year,agency) %>%
  summarise(CPA_120 = sum(CPA_120)) %>%
  group_by(agency)%>%
  mutate(CPA_120_cumul = cumsum(CPA_120))
d_CPA <- d_CPA[complete.cases(d_CPA),]
  
p_CPA <-
  ggplot(d_CPA,aes(x=fire_year,y=CPA_120_cumul,col=agency))+
  geom_line()+
  scale_color_manual(values=c("darkred","darkgreen","orange"))+
  labs(x = "fire year", y = "cumulative area >120 m in from patch edge (ha)")+
  scale_x_continuous(minor_breaks = seq(1984, 2015, 1),breaks = seq(1985, 2015, 5))+
  theme_bw()+
  theme(panel.grid.minor = element_line(size=0.5),
        panel.grid.major.x = element_line(color = "black"))

area_burned <- 
  group_by(d,agency) %>%
  summarise(area_burned=sum(firesize_ha))
d_CPA[d_CPA$fire_year==2015,"CPA_120_cumul"]/area_burned[c(1:3),"area_burned"]

png(file = paste0("./Figures/Fig5_",Sys.Date(),".png"),width=7.5,height=4,units="in",res=500)
p_CPA
dev.off()
  