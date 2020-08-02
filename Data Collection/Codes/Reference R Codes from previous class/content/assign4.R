#Econ 835 
#Assignment 3


# load vars package 
library(vars)
#load urca package
library(urca)


setwd("e:/835/prg_r") # Set my working directory
sink("assign3output.txt", append=TRUE, split=TRUE) #Save output from this script



# Question 1 
# data 1976.02 - 1996.06
data=read.table("usuk.txt",sep="",header=TRUE)


# 1 - get time series
dspot=diff(data[,"USUKS"])
fp=data[-1,"USUKF"]-data[-1,"USUKS"]
# full sample
vardata0=ts(cbind(dspot,fp),start=c(1976,3),frequency=12)


# 3 - find optimal number of lags
info.crit <- VARselect(vardata0,lag.max=4,type="const")
info.crit

# estimate VAR(1)
model0=VAR(vardata0,p=1,type="const")
summary(model0)




# get IRF
###Impact of shock to dspot
irf.dspot=irf(model0,impulse="dspot",n.ahead=24,ci=.9)
###Impact of shock to fp
irf.fp=irf(model0,impulse="fp",n.ahead=24,ci=.9)
###SVAR

amat=diag(2)
amat[2,1]=NA
amat
bmat=diag(2)
diag(bmat)=NA
bmat

svar0=SVAR(model0, estmethod = "direct", Amat = amat, Bmat = bmat,
hessian = TRUE, method = "BFGS")

irf.fp.svar=irf(svar0,impulse="fp",n.ahead=24,ci=.9)

irf.dspot.svar=irf(svar0,impulse="dspot",n.ahead=24,ci=.9)



#plot IRF (Comparing reduced form VAR IRF with structural VAR)
pdf(file = "graphsassign3.1.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(2,1))
## draw the plot
plot(irf.fp)
plot(irf.fp.svar)

## close the device to do the drawing
dev.off()



##########IMPOSING ORTHOGONALITY ASSUMPTION#############
amat1=diag(2)
amat1
bmat1=diag(2)
diag(bmat1)=NA
bmat1

svar1=SVAR(model0, estmethod = "direct", Amat = amat1, Bmat = bmat1,
hessian = TRUE, method = "BFGS")

irf.fp.svar1=irf(svar1,impulse="fp",n.ahead=24,ci=.9)



#plot IRF (Comparing reduced form VAR IRF with structural VAR)
pdf(file = "graphsassign3.2.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(2,1))
## draw the plot
plot(irf.fp)
plot(irf.fp.svar1)

## close the device to do the drawing
dev.off()





#########BQ Decomposition###########


data <-read.table("BQ.txt",sep="",header=TRUE)

# get time series
yt <- 100*log(data$GNP82)
dyt <- diff(yt)
unt <- data$USUNRATEE[-1]
nt <- length(dyt)


var.bq=ts(cbind(dyt,unt),start=c(1951,2),frequency=4)


# 3 - find optimal number of lags
info.bq <- VARselect(var.bq,lag.max=4,type="const")
info.bq

# estimate VAR
model.bq.rf=VAR(var.bq,p=2,type="const")
summary(model.bq.rf)


model.bq=BQ(model.bq.rf)
summary(model.bq)

# get IRF
irf.dy.bq=irf(model.bq,impulse="dyt",boot=FALSE,n.ahead=40)
irf.u.bq=irf(model.bq,impulse="unt",boot=FALSE,n.ahead=40)




# get supply shocks (dyt) 
# to get shocks to yt we need the cumulative sum 
supply <- cbind(cumsum(irf.dy.bq$irf$dyt[,1]),
							irf.dy.bq$irf$dyt[,2])

# get demand shocks (unt) 
# to get shocks to yt we need the cumulative sum 
# to match BQ graphs multiply by -1 
# (negative shock rahter than positive)
demand <- cbind(-1*cumsum(irf.u.bq$irf$unt[,1]),
							-1*irf.u.bq$irf$unt[,2])

#plot IRF (Comparing reduced form VAR IRF with structural VAR)
pdf(file = "graphsassign3.3.pdf")
## set up the new plotting device (pdf)
plot(demand[,1],type="l",col="black",lwd=2,ann=FALSE,
				xlim=c(0,24),ylim=c(-.6,1.4))
lines(demand[,2],col="blue",lwd=2)
abline(h=0)
legend(x="topright",c("Output Response", "Unemployment response"),
					col=c("black","blue"),lwd=2)


## close the device to do the drawing
dev.off()




# detrend data
zu <- cbind(matrix(1,nt,1),matrix(seq(1:nt)))


# function for OLS detrending
olsd <- function(yt,z){
	# construct ols detrended series
	bhat <- solve(t(z)%*%z)%*%t(z)%*%yt
	yd <- yt - z%*%bhat
	return(yd=yd)
	}

unt.1 <- olsd(unt,zu)

# set up data for VAR
var.bq1=ts(cbind(dyt,unt.1),start=c(1951,2),frequency=4)



info.bq.1 <- VARselect(var.bq1,lag.max=4,type="const")
info.bq.1

# estimate VAR
model.bq1.rf=VAR(var.bq,p=2,type="const")
summary(model.bq1.rf)


model.bq.1=BQ(model.bq1.rf)
summary(model.bq.1)

# get IRF
irf.dy.bq1=irf(model.bq.1,impulse="dyt",boot=FALSE,n.ahead=40)
irf.u.bq1=irf(model.bq.1,impulse="unt",boot=FALSE,n.ahead=40)




# get supply shocks (dyt) 
# to get shocks to yt we need the cumulative sum 
supply.1 <- cbind(cumsum(irf.dy.bq1$irf$dyt[,1]),
							irf.dy.bq1$irf$dyt[,2])

# get demand shocks (unt) 
# to get shocks to yt we need the cumulative sum 
# to match BQ graphs multiply by -1 
# (negative shock rahter than positive)
demand.1 <- cbind(-1*cumsum(irf.u.bq1$irf$unt[,1]),
							-1*irf.u.bq1$irf$unt[,2])


#plot IRF (Comparing reduced form VAR IRF with structural VAR)
pdf(file = "graphsassign3.4.pdf")
## set up the new plotting device (pdf)
plot(demand.1[,1],type="l",col="black",lwd=2,ann=FALSE,
				xlim=c(0,24),ylim=c(-.6,1.4))
lines(demand.1[,2],col="blue",lwd=2)
abline(h=0)
legend(x="topright",c("Output Response", "Unemployment response"),
					col=c("black","blue"),lwd=2)


## close the device to do the drawing
dev.off()


# compute FEVD
fevd.bq=fevd(model.bq,n.ahead=40)
fevd.bq1=fevd(model.bq.1,n.ahead=40)

#plot FEVD
pdf(file = "graphsassign3.5.pdf")
par(mfrow = c(2,1))

plot(fevd.bq)
plot(fevd.bq1)
## close the device to do the drawing
dev.off()

