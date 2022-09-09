rm(list=ls(all=TRUE))
ls()

# load dlm package 
library(dlm)
# load numDeriv package 
library(numDeriv)
# recession shading
library(tis)         

# Set seed
set.seed(123)

# Set directory
setwd("c:/835/prg_r") # Set my working directory

sink("assign5output.txt", append=TRUE, split=TRUE) #Save output from this script



(1)
# get data 
data <-read.table("tvp.txt",sep="",header=FALSE)
m1 <- ts(data$V2,start=c(1959,3),frequency=4)
dint <- ts(data$V3,start=c(1959,3),frequency=4)
inf <- ts(data$V4,start=c(1959,3),frequency=4)
surpl <- ts(data$V5,start=c(1959,3),frequency=4)
m1lag <- ts(data$V6,start=c(1959,3),frequency=4)

nobs <- length(m1)


X=cbind(dint,inf,surpl,m1lag)	
# TVP CAPM
# set parameter restrictions (only variances here)
parm_rest <- function(parm){
 	return( exp(parm) ) }

# set up SS model
ssm2 <- function(parm,x.mat){
	parm <- parm_rest(parm)
	return( dlmModReg(X=x.mat, dV=parm[1], 
dW=c(parm[2],parm[3],parm[4],parm[5],parm[6])) )
	}
# estimate parameters
fit2 <- dlmMLE(y=m1,parm=c(0.5,0.1,0.1,0.1,0.1,0.1),x.mat=X,build=ssm2,hessian=T)

# get estimates
se2 <- parm_rest(fit2$par)
se2

# get parameter estimates over time
# these are the smoothed state values
mod2 <- ssm2(fit2$par,X)
mod2f <- dlmFilter(m1,mod2)
mod2s <- dlmSmooth(mod2f)

# filtered estimates
beta0t <- ts(mod2f$m[-1,1],start=c(1959,3),frequency=4)
beta1t <- ts(mod2f$m[-1,2],start=c(1959,3),frequency=4)
beta2t <- ts(mod2f$m[-1,3],start=c(1959,3),frequency=4)
beta3t <- ts(mod2f$m[-1,4],start=c(1959,3),frequency=4)
beta4t <- ts(mod2f$m[-1,5],start=c(1959,3),frequency=4)

beta0t.1=ts(beta0t[11:length(beta0t)],start=c(1961,1),frequency=4)
beta1t.1=ts(beta1t[11:length(beta0t)],start=c(1961,1),frequency=4)
beta2t.1=ts(beta2t[11:length(beta0t)],start=c(1961,1),frequency=4)
beta3t.1=ts(beta3t[11:length(beta0t)],start=c(1961,1),frequency=4)
beta4t.1=ts(beta4t[11:length(beta0t)],start=c(1961,1),frequency=4)

#Saving the plot
pdf(file = "graphsassign4.5.pdf")

# plot parameters

par(mfrow = c(3,2))


plot(beta0t.1,plot.type='s',col=c("blue"),main="beta0t(filtered)")
plot(beta1t.1,plot.type='s',col=c("blue"),main="beta1t(filtered)")
plot(beta2t.1,plot.type='s',col=c("blue"),main="beta2t(filtered)")
plot(beta3t.1,plot.type='s',col=c("blue"),main="beta3t(filtered)")
plot(beta4t.1,plot.type='s',col=c("blue"),main="beta4t(filtered)")
## close the device to do the drawing
dev.off()


# smoothed estimates
beta0s <- ts(mod2s$s[-1,1],start=c(1959,3),frequency=4)
beta1s <- ts(mod2s$s[-1,2],start=c(1959,3),frequency=4)
beta2s <- ts(mod2s$s[-1,3],start=c(1959,3),frequency=4)
beta3s <- ts(mod2s$s[-1,4],start=c(1959,3),frequency=4)
beta4s <- ts(mod2s$s[-1,5],start=c(1959,3),frequency=4)


#Saving the plot
pdf(file = "graphsassign4.6.pdf")

# plot parameters

par(mfrow = c(3,2))


plot(beta0s,plot.type='s',col=c("blue"),main="beta0s(smoothed)")
plot(beta1s,plot.type='s',col=c("blue"),main="beta1s(smoothed)")
plot(beta2s,plot.type='s',col=c("blue"),main="beta2s(smoothed)")
plot(beta3s,plot.type='s',col=c("blue"),main="beta3s(smoothed)")
plot(beta4s,plot.type='s',col=c("blue"),main="beta4s(smoothed)")

## close the device to do the drawing
dev.off()


(2)

# Clark (1987,1989)

data <-read.table("rgdp_us.txt",sep="",header=FALSE)
yt <- ts(100*log(data$V2[1:251]),start=1948,frequency=4)
# US Unemployment rate
data <-read.table("ur_us.txt",sep="",header=FALSE)
ut <- ts(data$V2,start=c(1948,1),end=c(2010,3),frequency=4)

# get data 
yy <- cbind(yt,ut)
nobs <- length(yt)

# set parameter restrictions 
parm_rest <- function(parm){
 	parm[c(3,4,5,8)] <- exp(parm[c(3,4,5,8)])
 	return( parm ) 
	}

# set up SS model
ssm1 <- function(parm){
	parm <- parm_rest(parm)
	
	F.mat <- matrix(rep(0,10),nr=2)
	F.mat[1,c(1,3)] <- F.mat[2,5] <- 1
	F.mat[2,3:4] <- parm[6:7]
	V.mat <- diag(c(1e-07,parm[8]))
	
	G.mat <- matrix(rep(0,25),nr=5)
	G.mat[1,1:2] <- G.mat[2,2] <- G.mat[4,3] <- G.mat[5,5] <- 1
	G.mat[3,3:4] <- parm[1:2]
	W.mat <- diag(c(parm[4],0,parm[3],0,parm[5]))

	m0.mat <- matrix(rep(0,5),nr=5)
	C0.mat <- diag(5)*10^7
	
	return( dlm(FF=F.mat,V=V.mat,GG=G.mat,W=W.mat,
				m0=m0.mat,C0=C0.mat) )
	}

# estimate parameters
parm.start <- c(1.53,-.57,-0.8439,-1.0788,-1,-.34,-.16,-3.5065)
fit1 <- dlmMLE(y=yy,parm=parm.start,build=ssm1,hessian=T)
mod1 <- ssm1(fit1$par)

# filter and smooth
mod1f <- dlmFilter(yy,mod1); mod1s <- dlmSmooth(mod1f)

# get estimates for ARMA(2,0) part
coef <- parm_rest(fit1$par)
coef

# get parameter estimates 
drift <- mod1f$m[nobs+1,2]
covar <- dlmSvd2var(mod1f$U.C[[nobs+1]],mod1f$D.C[nobs+1,])
coef.se <- sqrt(covar[2,2])
drift; coef.se

# smoothed values
yt.trend <- ts(mod1s$s[-1,1],start=1948,frequency=4)
yt.cycle <- ts(mod1s$s[-1,3],start=1948,frequency=4)
ut.trend <- ts(mod1s$s[-1,5],start=1948,frequency=4)
ut.cycle <- ut - ut.trend

#Saving the plot
pdf(file = "graphsassign4.7.pdf")


# yt: plot filtered states (trend and cycle)
plot(yt.cycle,ylim=c(-10,7),xlim=c(1948,2011))
nberShade()
lines(yt.cycle)
abline(h=0)
dev.off()
#Saving the plot
pdf(file = "graphsassign4.8.pdf")


plot(cbind(yt,yt.trend),ylim=c(740,960),xlim=c(1948,2011),
			plot.type='s',col=c("black","blue"))
nberShade()
lines(yt)
lines(yt.trend)
dev.off()

#Saving the plot
pdf(file = "graphsassign4.9.pdf")
# ut: plot smoothed state (trend and cycle)
plot(ut.cycle,ylim=c(-4,5),xlim=c(1948,2011))
nberShade()
lines(ut.cycle)
abline(h=0)
dev.off()

#Saving the plot
pdf(file = "graphsassign4.10.pdf")

plot(cbind(ut,ut.trend),ylim=c(0,12),xlim=c(1948,2011),
			plot.type='s',col=c("black","blue"))

dev.off()

