
#packages
library(INLA)
library(reshape2)
library(rgeos)
library(fields)

  
#preparation - mesh construction - use the loc.domain argument

mesh <- inla.mesh.2d(loc.domain = dat2[,c(9,10)],max.edge=c(10000,25000),cutoff=8000, offset = c(50000,100000))

#plot the mesh to see what it looks like
plot(mesh)

##set the spde representation to be the mesh just created
spde <- inla.spde2.matern(mesh)

#make A matrix for structured data - should this be pulling the x and y coordinates for the location?
EA_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(EA_data[,9:10]))

#make A matrix for unstructured data
CS_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(CS_data[,5:6]))


# Joint model

# One spatial field
# Uses Simpson approach for PP data
# Binomial model for PA data
# Using cloglog


# create integration stack

#make dual mesh
dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2])
tiles <- deldir::tile.list(dd)

#make domain into spatial polygon
domainSP <- SpatialPolygons(list(Polygons(
  list(Polygon(ENG)), '0')))

#ensure there are no self intersections using a 0 buffer
domSP <- gBuffer(domainSP, byid=TRUE, width=0)
  
  
#intersection between domain and dual mesh
poly.gpc <- as(domSP@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")

# w now contains area of voronoi polygons
w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc)))

#check some have 0 weight
table(w>0)

#################################################
#################################################


nv <- mesh$n
n <- nrow(CS_data)


#change data to include 0s for nodes and 1s for presences
y.pp <- rep(0:1, c(nv, n))

#add expectation vector (area for integration points/nodes and 0 for presences)
e.pp <- c(w, rep(0, n))

#diagonal matrix for integration point A matrix
imat <- Diagonal(nv, rep(1, nv))

A.pp <- rBind(imat, CS_data_A)

#get covariate for integration points


mesh.sp.prox <- paste(((1000 * ceiling(mesh$loc[,1]/1000))-500),((1000 * ceiling(mesh$loc[,2]/1000))-500),sep="_")
msh.mt <- match(mesh.sp.prox,grid1km$sp.prox)
sources.mesh = grid1km$sources[msh.mt]
nodes.mesh = grid1km$nodes[msh.mt]
urban.mesh = grid1km$Urban[msh.mt]
sources_node_ratio.mesh = grid1km$sources_node_ratio[msh.mt]

#data stack for point based cit sci data

stk_CS_data <- inla.stack(data=list(y=cbind(y.pp, NA), e = e.pp),
                      effects=list(list(data.frame(interceptB=rep(1,nv+n)),sources=c(sources.mesh,CS_data$sources),nodes=c(nodes.mesh,CS_data$nodes),Urban=c(urban.mesh,CS_data$Urban),sources_node_ratio=c(sources_node_ratio.mesh,CS_data$sources_node_ratio)), list(uns_field=1:spde$n.spde, bias_field = 1:spde$n.spde)),
                      A=list(1,A.pp),
                      tag="CS_data")	

					  
#data stack for structured EA data

stk_EA_data <- inla.stack(data=list(y=cbind(NA, EA_data$Pres), Ntrials = rep(1, nrow(EA_data))),
                      effects=list(list(data.frame(interceptA=rep(1,length(EA_data$FULL_EASTING))),sources=EA_data$sources,nodes=EA_data$nodes,Urban=EA_data$Urban,sources_node_ratio=EA_data$sources_node_ratio), list(str_field=1:spde$n.spde)),
                      A=list(1,EA_data_A),
                      tag="EA_data")
					  
					  
##NOTE: doesn't use the copy function initially

stk <- inla.stack(stk_CS_data, stk_EA_data)


###################################
###################################

# Create a grid to use for predictions - IGNORED IN THE END AND MODEL FORMULA USED
  
ys <- cbind(rep(NA, nrow(grid1km)), rep(NA, nrow(grid1km)))

A.pred <- inla.spde.make.A(mesh, loc=as.matrix(grid1km[,1:2]))

stack.pred <- inla.stack(data=list(y=ys),
				effects = list(list(data.frame(interceptA=rep(1,nrow(grid1km)))), sources=grid1km$sources,nodes=grid1km$nodes,Urban=grid1km$Urban, list(uns_field=1:spde$n.spde)),
           A=list(1,1,1,1, A.pred),
           tag='pred')
  

#join.stack <- inla.stack(stk,stack.pred)

join.stack <- stk

########################################
########################################
########################################

formulaJ = y ~  interceptA + interceptB + sources + nodes + Urban + sources_node_ratio + f(uns_field, model = spde) + f(str_field, copy = "uns_field", fixed = TRUE) + f(bias_field,model=spde) -1


result <- inla(formulaJ,family=c("poisson", "binomial"),
               data=inla.stack.data(join.stack),
               control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE),
               control.family = list(list(link = "log"), 
                                     list(link = "cloglog")),
               E = inla.stack.data(join.stack)$e,
               Ntrials = inla.stack.data(join.stack)$Ntrials,
               control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
)


########################################
########################################
########################################


##project the mesh onto the initial simulated grid 100x100 cells in dimension
proj1 <- inla.mesh.projector(mesh,ylim=c(0,max(ENG[,2])),xlim=c(1,max(ENG[,1])),dims=c(max(ENG[,1]/1000),max(ENG[,2])/1000))

proj.df <- inla.mesh.projector(mesh,loc=as.matrix(grid1km[,1:2]))


##pull out the mean of the random field for the CS model
bmean <- inla.mesh.project(proj1, result$summary.random$bias_field$mean)
smean <- inla.mesh.project(proj1, result$summary.random$str_field$mean)
umean <- inla.mesh.project(proj1, result$summary.random$uns_field$mean)

##plot the standard deviation of random field
bsd <- inla.mesh.project(proj1, result$summary.random$bias_field$sd)
ssd <- inla.mesh.project(proj1, result$summary.random$str_field$sd)
usd <- inla.mesh.project(proj1, result$summary.random$uns_field$sd)


par(mfrow=c(1,2))
image.plot(seq(1,max(ENG[,1]),by=1000),seq(1,max(ENG[,2]),by=1000),bmean)
lines(ENG,col="white",lwd=3)
image.plot(seq(1,max(ENG[,1]),by=1000),seq(1,max(ENG[,2]),by=1000),bsd)
lines(ENG,col="white",lwd=3)


res=result$summary.fixed




grid1km$rnfld <- inla.mesh.project(proj.df, result$summary.random$str_field$mean)
grid1km$rnfld.sd <- inla.mesh.project(proj.df, result$summary.random$str_field$sd)


# Use the model to produce estimates of the number of cases using the spatial field and ignoring it. 
grid1km$Est = (res[1,1] + apply(t(res[3:6,1]*t(grid1km[,c(4,5,8,6)])),1,sum) + grid1km$rnfld)
grid1km$pred.prob=(1-exp(-exp(grid1km$Est)))

grid1km$ENG <- inout(grid1km[,1:2],ENG)

colr=viridis(10)
#plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[ceiling((grid1km$Est[grid1km$ENG]*10)+104)],cex=0.15,asp=1)
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[as.integer(as.factor(cut(grid1km$Est[grid1km$ENG],c(-200,-8:0,200))))],cex=0.15,asp=1)

library(viridis)
colr=viridis(10)

library(splancs)
GB2 <- GB[complete.cases(GB),]
in.GB <- inout(occ[,5:6],GB2)
occGB <- occ[in.GB,]

png("Figure 3 - 3 panel version.png", height = 1600, width = 2800, res = 300)
par(mfrow=c(1,3), xpd = TRUE)
plot(dat1yr[,9:10],asp=1,xlab="",ylab="", xlim = c(0, 700000), ylim = c(0, 1000000), frame= FALSE, axes = FALSE, pch = 20, col = "grey")
text(600000,1200000, "a)", cex = 2, font = 2)
points(dat2yr[,9:10],pch=20,col="black")
lines(GB)
plot(occGB[,5:6],asp=1,xlab="",ylab="", xlim = c(0, 700000), ylim = c(0, 1000000), frame= FALSE, axes = FALSE, pch = 20, col = "blue")
text(600000,1200000, "b)", cex = 2, font = 2)
lines(GB)
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[as.integer(as.factor(cut(grid1km$Est[grid1km$ENG],c(-200,-8:0,200))))],cex=0.15,asp=1,xaxt="n",yaxt="n", xlim = c(0, 700000), ylim = c(0, 1000000), ylab = "", xlab = "", frame = FALSE, axes = FALSE)
legend(0, 1200000,title=as.expression(bquote(bold("Occurrence probability\n (on cloglog scale)"))),bty="n",legend=rev(c("< -8","-8:-7","-7:-6","-6:-5","-5:-4","-4:-3","-3:-2","-2:-1","-1:0","> 0")),pch=15,col=rev(colr),cex=1.2,pt.cex=1.3)
text(600000,1200000, "c)", cex = 2, font = 2)
lines(ENG)
#legend("bottomright",title="Integrated model - mean estimate",bty="n",legend="")
dev.off()

############################

pboot <- matrix(ncol=100,nrow=length(grid1km[,1]))

for(k in 1:100){

	int.boot <- rnorm(1,res[1,1],res[1,2])
	coeff.boot <- c(rnorm(1,res[3,1],res[3,2]),rnorm(1,res[4,1],res[4,2]),rnorm(1,res[5,1],res[5,2]),rnorm(1,res[6,1],res[6,2]))
	rf.boot <- rnorm(length(grid1km[,1]),grid1km$rnfld,grid1km$rnfld.sd)


	# Use the model to produce estimates of the number of cases using the spatial field and ignoring it. 
	bt.linpred = (int.boot + apply(t(coeff.boot*t(grid1km[,c(4,5,8,6)])),1,sum) + rf.boot)
	pboot[,k]=(1-exp(-exp(bt.linpred)))

}


grid1km$err <- apply(pboot,1,sd)
grid1km$coef.var <- (grid1km$pred.prob/grid1km$err)

colr=grey(seq(0.05,0.95,len=50))
colr=c(colr,rep(colr[50],5000))
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[ceiling(((grid1km$coef.var[grid1km$ENG])*100))],cex=0.15,asp=1)






