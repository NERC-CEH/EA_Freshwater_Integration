stk_EA_data <- inla.stack(data=list(y=EA_data$Pres, Ntrials = rep(1, nrow(EA_data))),
                      effects=list(list(data.frame(interceptA=rep(1,length(EA_data$FULL_EASTING))),sources=EA_data$sources,nodes=EA_data$nodes,Urban=EA_data$Urban,sources_node_ratio=EA_data$sources_node_ratio), list(str_field=1:spde$n.spde)),
                      A=list(1,EA_data_A),
                      tag="EA_data")

join.stack <- stk_EA_data

########################################
########################################
########################################

formulaJ = y ~  interceptA + sources + nodes + Urban + sources_node_ratio + f(str_field,model=spde) -1


result.EA <- inla(formulaJ,family=c("binomial"),
               data=inla.stack.data(join.stack),
               control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE),
               control.family = list(link = "logit"),
               Ntrials = inla.stack.data(join.stack)$Ntrials,
               control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
)

res.EA=result.EA$summary.fixed

########################################
########################################
########################################

grid1km$rnfld.EA <- inla.mesh.project(proj.df, result.EA$summary.random$str_field$mean)
grid1km$rnfld.sd.EA <- inla.mesh.project(proj.df, result.EA$summary.random$str_field$sd)


# Use the model to produce estimates of the number of cases using the spatial field and ignoring it. 
grid1km$Est.EA = (res.EA[1,1] + apply(t(res.EA[2:5,1]*t(grid1km[,c(4,5,8,6)])),1,sum) + grid1km$rnfld.EA)
grid1km$pred.prob.EA=(1-exp(-exp(grid1km$Est.EA)))


colr=tim.colors(100)
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[ceiling((grid1km$Est.EA[grid1km$ENG]*10)+100)],cex=0.15)


############################

pboot.EA <- matrix(ncol=100,nrow=length(grid1km[,1]))

for(k in 1:100){

	int.boot <- rnorm(1,res.EA[1,1],res.EA[1,2])
	coeff.boot <- c(rnorm(1,res.EA[2,1],res.EA[2,2]),rnorm(1,res.EA[3,1],res.EA[3,2]),rnorm(1,res.EA[4,1],res.EA[4,2]),rnorm(1,res.EA[5,1],res.EA[5,2]))
	rf.boot <- rnorm(length(grid1km[,1]),grid1km$rnfld.EA,grid1km$rnfld.sd.EA)


	# Use the model to produce estimates of the number of cases using the spatial field and ignoring it. 
	bt.linpred = (int.boot + apply(t(coeff.boot*t(grid1km[,c(4,5,8,6)])),1,sum) + rf.boot)
	pboot.EA[,k]=(1-exp(-exp(bt.linpred)))

}


grid1km$err.EA <- apply(pboot.EA,1,sd)
summary(grid1km$err.EA)
grid1km$coef.var.EA <- (grid1km$pred.prob.EA/grid1km$err.EA)

colr=rev(grey(seq(0.05,0.95,len=12)))
colr=c(colr,rep(colr[12],5000))
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[ceiling(((grid1km$err.EA[grid1km$ENG])*100))],cex=0.15)




######################
######################
colr=tim.colors(10)
#plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[ceiling((grid1km$Est[grid1km$ENG]*10)+104)],cex=0.15,asp=1)
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[as.integer(as.factor(cut(grid1km$Est[grid1km$ENG],c(-200,-8:0,200))))],cex=0.15,asp=1)

par(mfrow=c(2,2),mai=c(0.1,0.2,0.1,0.1))

colr=tim.colors(10)
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[as.integer(as.factor(cut(grid1km$Est[grid1km$ENG],c(-200,-8:0,200))))],cex=0.15,asp=1,xaxt="n",yaxt="n")
legend("topleft",title="Occurrence probability (on cloglog scale)",bty="n",legend=c("< -8","-8:-7","-7:-6","-6:-5","-5:-4","-4:-3","-3:-2","-2:-1","-1:0","> 0"),pch=15,col=colr,cex=0.7,pt.cex=1.3)
legend("bottomright",title="Integrated model - mean estimate",bty="n",legend="")

colr=rev(grey(seq(0.05,0.95,len=12)))
colr=c(colr,rep(colr[12],5000))
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[ceiling(((grid1km$err[grid1km$ENG])*100))],cex=0.15,asp=1,xaxt="n",yaxt="n")
lines(ENG)
legend("bottomright",title="Integrated model - uncertainty",bty="n",legend="")
legend("topleft",title="Standard Error of prediction",bty="n",legend=c("0:0.01","0.01:0.02","0.02:0.03","0.03:0.04","0.04:0.05","0.05:0.06","0.06:0.07","0.07:0.08","0.08:0.09","0.09:0.1","> 0.1"),pch=15,col=colr[1:11],cex=0.7,pt.cex=1.3)



colr=tim.colors(10)
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[as.integer(as.factor(cut(grid1km$Est.EA[grid1km$ENG],c(-200,-8:0,200))))],cex=0.15,asp=1,xaxt="n",yaxt="n")
legend("topleft",title="Occurrence probability (on cloglog scale)",bty="n",legend=c("< -8","-8:-7","-7:-6","-6:-5","-5:-4","-4:-3","-3:-2","-2:-1","-1:0","> 0"),pch=15,col=colr,cex=0.7,pt.cex=1.3)
legend("bottomright",title="Agency only model - mean estimate",bty="n",legend="")


colr=rev(grey(seq(0.05,0.95,len=12)))
colr=c(colr,rep(colr[12],5000))
plot(grid1km[grid1km$ENG,1:2],pch=15,col=colr[ceiling(((grid1km$err.EA[grid1km$ENG])*100))],cex=0.15,asp=1,xaxt="n",yaxt="n")
lines(ENG)
legend("bottomright",title="Agency only model - uncertainty",bty="n",legend="")
legend("topleft",title="Standard Error of prediction",bty="n",legend=c("0:0.01","0.01:0.02","0.02:0.03","0.03:0.04","0.04:0.05","0.05:0.06","0.06:0.07","0.07:0.08","0.08:0.09","0.09:0.1","> 0.1"),pch=15,col=colr[1:11],cex=0.7,pt.cex=1.3)


