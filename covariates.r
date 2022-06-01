###################################
## 
##
##
##
###################################

#load libraries required
library(shapefiles)
library(reshape2)
library(caret)
library(gdata)
library(rgeos)
library(sp)
library(rgdal)
library(sf)
library(raster)


## create a 1km grid over the domain and store both as sf object and simple data frame

r <- raster(extent(matrix( c(130000, 660000, 0, 660000), byrow=T,nrow=2)), nrow=660, ncol=530, 
            crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs") 
r[] <- 1:ncell(r)
r2 = as(r, "SpatialPolygonsDataFrame")
r3 = st_as_sf(r2,crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs")

grid1km <- data.frame(st_coordinates(st_centroid(r3)))
grid1km$Layer <- 1:length(grid1km[,1])



################################
###### River network data ########
################################


## read in river network covariate data - nodes and sources

cpm = read_sf("C://TempData//EA//IRN_Sources.shp")
cpm_points = st_centroid(cpm) 
cpm_points = st_transform(cpm_points, crs="+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs")

cpm2 = read_sf("C://TempData//EA//IRN_Nodes.shp")
cpm2_points = st_centroid(cpm2) 
cpm2_points = st_transform(cpm2_points, crs="+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs")


## intesect the 1km grid and river network data to provide spatial covariate layers

nc_point_in_poly <- data.frame(st_join(cpm_points, r3, join = st_within))
nc_point_in_poly2 <- data.frame(st_join(cpm2_points, r3, join = st_within))

##count the sources and nodes for each pixel in the 1km raster - this defines the spatial covariates to use
sources1km <- aggregate(nc_point_in_poly$OBJECTID~nc_point_in_poly$layer,FUN=length)
nodes1km <- aggregate(nc_point_in_poly2$OBJECTID~nc_point_in_poly2$layer,FUN=length)

#match the corresponding pixels from the 1km grid and add to the data frame created earlier to give full raster covariates
idx <- match(sources1km[,1],grid1km$Layer)
grid1km$sources=0
grid1km$sources[idx]=sources1km[,2]

idx.n <- match(nodes1km[,1],grid1km$Layer)
grid1km$nodes=0
grid1km$nodes[idx.n]=nodes1km[,2]


grid1km$sources_node_ratio <- (grid1km$sources/grid1km$nodes)
grid1km$sources_node_ratio[grid1km$nodes==0]=0

#create a simple spatial proxy for use matching up later on
grid1km$sp.prox <- paste(grid1km$X,grid1km$Y,sep="_")



################################
###### Percentage Urban Area ########
################################

## read in the percent target class data from the land cover map. 
dat <- read.table("C://TempData//EA//LCM2007_GB_1K_PC_TargetClass_22.asc",skip=6)


##set up a data frame of spatial locations according to LCM ascii file
east=seq(0,by=1000,len=700)
north=seq(0,by=1000,len=1300)
lc_grid=expand.grid(East=east,North=north)

## create third column representing entries in ascii grid. This isnt very efficient. 
v=c()
for(i in 1300:1){
	v=c(v,(dat[i,]))
} 
lc_grid$Val=unlist(v)

#plot the resulting data frame values to ensure correct
plot(lc_grid[,1:2],pch=15,cex=0.1,col=grey(1-(lc_grid$Val/100)),asp=1)

#create a simple spatial proxy for use matching up later on
lc_grid$sp.prox <- paste(lc_grid$East+500,lc_grid$North+500,sep="_")


#match the corresponding pixels from the 1km grid and add to the data frame created earlier to give full raster covariate for % urban
idx2 <- match(grid1km$sp.prox,lc_grid$sp.prox)
grid1km$Urban <- lc_grid$Val[idx2]


##############################
##############################

### match and append covariate data to the EA data and citizen science data sets previously created
EA.cov.id <- match(EA_data$sp.prox,grid1km$sp.prox)
EA_data=cbind(EA_data,grid1km[EA.cov.id,c(4,5,6,8)])

CS.cov.id <- match(CS_data$sp.prox,grid1km$sp.prox)
CS_data=cbind(CS_data,grid1km[CS.cov.id,c(4,5,6,8)])


##############################
##############################

grid1km$ENG <- inout(grid1km[,1:2],ENG)

par(mfrow=c(2,2),mai=c(0,0,0,0))

library(fields)
colr=tim.colors(20)
colr=c(colr,rep(colr[20],40))
plot(grid1km[grid1km$ENG,1:2],pch=15,cex=0.15,col=colr[grid1km$nodes[grid1km$ENG]+1],asp=1,xaxt="n",yaxt="n",xlab="",ylab="")
legend("topleft",col=colr[1:20],pch=15,pt.cex=2,legend=c(0:18,">18"),bty="n",title="Nodes per km")

 
colr=tim.colors(8)
colr=c(colr,rep(colr[8],40))
plot(grid1km[grid1km$ENG,1:2],pch=15,cex=0.15,col=colr[grid1km$sources[grid1km$ENG]+1],asp=1,xaxt="n",yaxt="n",xlab="",ylab="")
legend("topleft",col=colr[1:8],pch=15,pt.cex=2,legend=c(0:6,">6"),bty="n",title="Sources per km")
  
 
colr=tim.colors(11)
plot(grid1km[grid1km$ENG,1:2],pch=15,cex=0.15,col=colr[((grid1km$sources_node_ratio[grid1km$ENG])*10)+1],asp=1,xaxt="n",yaxt="n",xlab="",ylab="")
legend("topleft",col=colr,pch=15,pt.cex=2,legend=seq(0,1,by=0.1),bty="n",title="Sources to Nodes ratio")

colr=tim.colors(100)
plot(grid1km[grid1km$ENG,1:2],pch=15,cex=0.15,col=colr[((grid1km$Urban[grid1km$ENG]))+1],asp=1,xaxt="n",yaxt="n",xlab="",ylab="")
legend("topleft",col=colr[seq(1,100,by=10)],pch=15,pt.cex=2,legend=seq(0,99,by=10),bty="n",title="% Urban")






