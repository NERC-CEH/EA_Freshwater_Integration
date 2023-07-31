########################################################################
###
###    EA Freshwater integrated species distribution modelling    ######
###
###		
###
########################################################################

#setwd("C://TempData//EA//")

### Read in the agency data as downloaded from: https://environment.data.gov.uk/ecology-fish/downloads/
tx <- read.csv("INV_OPEN_DATA_TAXA.csv")
sp_id <- read.csv("TAXA_INFO.csv")
sites <- read.csv("INV_OPEN_DATA_SITE.csv")
metrics <- read.csv("INV_OPEN_DATA_METRICS.csv")


### remove all columns from data sets that are not required to keep data tidy 
metrics_red <- metrics[,c(1,2,11,12)]
sites_red <- sites[,c(1,2,5,6,12,13)]
tx_red <- tx[,c(4,7,8,9)]


#combine the metrics and the site data into merged data frame
dat=cbind(metrics_red,sites_red[match(metrics_red[,1],sites_red[,4]),])


#append the taxonomic data to this merged data frame for relevant site and sample
dat1=cbind(dat,tx_red[match(dat[,3],tx_red[,1]),])
#create seperate year column from the data string
dat1$Year=strptime(dat1$DATE_OF_ANALYSIS,"%d/%m/%Y")$year+1900


#filter the data to observations of water beetles of agabus genus
dat2=dat1[is.element(dat1$TAXON_LIST_ITEM_KEY,sp_id$TAXON_LIST_ITEM_KEY[sp_id$PARENT_TAXON_NAME=="Agabus"]),]



##############################################

##############################################


### read in citizen science (CS) recording scheme data on water beetles, downloaded via NBN
occ=read.table("occurrence.txt",header=T,fill=T,sep="\t")


#filter the data to records of agabus genus since 2000
occ=occ[occ$year>2000 & occ$genus=="Agabus",c(103,107,133:134)]
occ=occ[-which(is.na(occ$year)),]


#ensure that the spatial data has been read as a numeric rather than as a factor
occ[,3]=as.numeric(as.character(occ[,3]))
occ[,4]=as.numeric(as.character(occ[,4]))


#convert lat lon coordinates to BNG eastings and northing to be compatible with agency records' site descriptions
library(rgdal)
cord.dec = SpatialPoints(cbind(occ[,4], occ[,3]), proj4string = CRS("+proj=longlat"))
cord.BNG <- spTransform(cord.dec, CRS("+init=epsg:27700"))

#store the BNG coordinates in the same data frame
occ$EASTING=data.frame(cord.BNG)[,1]
occ$NORTHING=data.frame(cord.BNG)[,2]


#read in GB coast file to provide boundary polygon
GB <- read.table("GBoutline.txt",header=T)
#specify England polygon only
ENG <- GB[1:1662,]


#filter the CS occurrence data to that only from within England
library(splancs)
in.ENG <- inout(occ[,5:6],ENG)
occE <- occ[in.ENG,]



#############################

#############################



##select the agency data from 2001:2003 matching the CS data 
dat1yr <- dat1[dat1$Year>2000 & dat1$Year<2004,]
dat2yr <- dat2[dat2$Year>2000 & dat2$Year<2004,]

##create a presence/absence column for the EA data from the occurrence records 
dat1yr <- dat1yr[-which(duplicated(dat1yr[,1])),]
dat1yr$Pres <- 0
dat1yr$Pres[is.element(dat1yr$ï..SITE_ID,dat2yr$ï..SITE_ID)]=1


#rename final data frames for use later on 
EA_data=dat1yr
CS_data=occE

#create a spatial proxy column that will enable simple lookups later on
EA_data$sp.prox <- paste(((1000 * ceiling(EA_data[,9]/1000))-500),((1000 * ceiling(EA_data[,10]/1000))-500),sep="_")
CS_data$sp.prox <- paste(((1000 * ceiling(CS_data[,5]/1000))-500),((1000 * ceiling(CS_data[,6]/1000))-500),sep="_")



#############################

#############################



## simple plot of data to check all looks ok and there are no spurious records or issues 
par(mfrow=c(1,2))
plot(dat1yr[,9:10],main="EA Monitoring",asp=1,xlab="Easting",ylab="Northing")
points(dat2yr[,9:10],pch=20,col="red")
lines(ENG)
plot(occE[,5:6],main="Scheme Data",asp=1,xlab="Easting",ylab="Northing")
lines(ENG)


