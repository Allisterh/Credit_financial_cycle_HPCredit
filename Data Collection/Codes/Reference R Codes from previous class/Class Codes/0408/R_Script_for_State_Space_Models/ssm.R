# load dlm package 
library(dlm)
library(tis)
# Set directory
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0408/R_Script_for_State_Space_Models") 

# Set seed
set.seed(1234)


##################################################################
# SS Model
# 
# y(t) = FF a(t) + v(t)		v(t)~N(0,V(t))
# a(t) = GG a(t-1) + w(t)	w(t)~N(0,W(t))
#							a(0)~N(m0,C0)
# 
##################################################################


##################################################################
# Example 1: AR(1) 
##################################################################

# simulate AR(1) process
# phi = .8, sig2 = .25
nobs=250
yt=arima.sim(n=nobs,list(ar=.8,ma=0),sd=.5)
plot(yt)
# estimate AR(1) for comparison
model10=arima(yt,order=c(1,0,0),method="ML",include.mean=FALSE)
model10

# write a function
# set parameter restrictions 
parm_rest=function(parm){
 	return( c(exp(parm[1])/(1+exp(parm[1])),exp(parm[2])) ) 
	}

#exponential function always return positive values
# set up SS model
ssm1=function(parm){
	parm=parm_rest(parm)
	return( dlm(FF=1,V=0,GG=parm[1],W=parm[2],
				m0=0,C0=solve(1-parm[1]^2)*parm[2]) )
}

##FF and GG comments on notes, FF in 
# estimate parameters
        fit1=dlmMLE(y=yt,parm=c(0,1),build=ssm1,hessian=T)

# get estimates 
coef=parm_rest(fit1$par)
coef
# get standard errors using delta method
dg1=exp(fit1$par[1])/(1+exp(fit1$par[1]))^2
dg2=exp(fit1$par[2])
dg=diag(c(dg1,dg2))
var=dg%*%solve(fit1$hessian)%*%dg
# print results
coef; sqrt(diag(var))

# print SS model
ssm1(fit1$par)


# AR(1) again

# a simpler option is to use DLM package functions
# set up SS model
ssm2 =function(parm){
	parm=parm_rest(parm)
	dlm=dlmModARMA(ar=parm[1], ma=NULL, sigma2=parm[2])
	dlm$C0=solve(1-parm[1]^2)*parm[2] 
	return(dlm)
	}
# estimate parameters
fit2=dlmMLE(y=yt,parm=c(0,1),build=ssm2,hessian=T)
fit2

# get estimates 
coef=parm_rest(fit2$par)
# get standard errors using delta method
dg1=exp(fit2$par[1])/(1+exp(fit2$par[1]))^2
dg2=exp(fit2$par[2])
dg=diag(c(dg1,dg2))
var=dg%*%solve(fit2$hessian)%*%dg
# print results
coef; sqrt(diag(var))



##################################################################
# Example 2: TV-CAPM
# from Analysis of Financial Time Series (Ruey Tsay - 2010)
##################################################################

# get data (GM and SP500 excess returns)
# sample: 1990.1 - 2008.12
data=read.table("capm.txt",sep="",header=TRUE)
gm=ts(100*data$gm,start=c(1990,1),frequency=12)
sp=ts(100*data$sp,start=c(1990,1),frequency=12)
nobs=length(gm)

# fit least squares model
fit0=lm(gm~sp)
summary(fit0)

	
# TVP CAPM
# set parameter restrictions (only variances here)
parm_rest=function(parm){
 	return( exp(parm) ) }

# set up SS model
ssm2=function(parm,x.mat){
	parm=parm_rest(parm)
	return( dlmModReg(X=x.mat, dV=parm[1], dW=c(parm[2],parm[3])) )
}
  #impose restriction on the regression model

#X:the design matrix

#dV:variance of the observation noise.
#dW: diagonal elements of the variance matrix of the system noise.

# estimate parameters
fit2=dlmMLE(y=gm,parm=c(1,1,1),x.mat=sp,build=ssm2,hessian=T)

# get estimates
se2=parm_rest(fit2$par)
se2
  #estimate of epsilon_t, e1t, e2t in GM CAPM model

# get parameter estimates over time
# these are the smoothed state values
mod2=ssm2(fit2$par,sp) #unconstrainted para
mod2f=dlmFilter(gm,mod2) #constrained para
mod2s=dlmSmooth(mod2f)  #smoothed 

# filtered states
alphat=ts(mod2f$m[-1,1],start=c(1990,1),frequency=12) #intercept , -1 : throw away first row of 0s
betat=ts(mod2f$m[-1,2],start=c(1990,1),frequency=12) #slopes
# Smoothed states
betas=ts(mod2s$s[-1,2],start=c(1990,1),frequency=12)

# plot parameters
plot(alphat,plot.type='s',col=c("blue"),
		main="TVP-CAPM - Intercept (filtered estimate)")

plot(betat,plot.type='s',col=c("blue"),
		main="TVP-CAPM - Slope (filtered estimate)")
	
plot(betas,plot.type='s',col=c("blue"),
     main="TVP-CAPM - Slope (smoothed estimate)")
