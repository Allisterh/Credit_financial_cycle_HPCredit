rm(list=ls(all=TRUE))
ls()

# load urca package 
library(urca)
library(vars)
library(MASS)
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0506") # Set my working directory


# Set seed
set.seed(1234)


# Spurious regression
e1 <- rnorm(250)
e2 <- rnorm(250)
y1 <- cumsum(e1)
y2 <- cumsum(e2)

# get regression results with variables in levels
summary(lm(y1~y2))

# get regression results with variables in diff
summary(lm(diff(y1)~diff(y2)))


# Example bivariate cointegrated system 
e1 <- rnorm(250,mean=0,sd=.5)
e2 <- rnorm(250,mean=0,sd=.5)
u.ar1 <- arima.sim(model=list(ar=.75),250,innov=e1) #shock is e1
y2 <- ts(cumsum(e2))
y1 <- ts(y2 + u.ar1)
  #setting up cointegration

# plot time series
layout(1)
plot(y1,col="blue")
lines(y2)
# plot cointegrating residual
plot(u.ar1,col="blue")
abline(h=0)


# Example Zivot (2000)
# US-UK data 1976.02 - 1996.06
data <-read.table("usuk.txt",sep="",header=TRUE)

# get time series
st <- ts(data[,"USUKS"],start=c(1976,2),frequency=12)
ft <- ts(data[,"USUKF"],start=c(1976,2),frequency=12)

# plot time series
plot(st,col="blue")
lines(ft)


# test for cointegration with pre-specified cointegrating vector
# get interest rate differential
ut <- ts(st-ft,start=c(1976,2),frequency=12)

# plot cointegrating vector , difference in those two series, spot and forward rates
plot(ut,col="blue")
abline(h=0)
# test for unit root in cointegrating residual
ur.test=ur.pp(ut,type = c("Z-tau"), model = c("constant"))
summary(ur.test)

# test for cointegration with estimated cointegrating vector
coint.eq <- lm(st~ft)
summary(coint.eq) #this shows that the integration vector is (1,-1)
ut.est <- coint.eq$resid

# test for unit root in cointegrating residual 
ut.ur=ur.pp(ut.est,type = c("Z-tau"), model = c("constant"))
summary(ut.ur)

#Phillips-Ouliaris test
po.test=ca.po(cbind(st,ft),demean= "none", type ="Pu")
summary(po.test)
#reject this test of unitroot-> cointegration
#two types Pu and Pz to choose here

#VECM 
  #to find lags : do VAR in levels

#########VECM ANALYSIS####

vardata0=ts(cbind(ft,st),start=c(1976,2),frequency=12)
info.crit=VARselect(vardata0,lag.max=8,type="const")
info.crit #-> this shows we only need one lag<- AIC value


dst=na.omit(diff(st))
dft=na.omit(diff(ft))
ect.1=ut.est[2:(length(ut.est)-1)]
ect.pre=ut[2:(length(ut)-1)]


dft.t=dft[2:(length(dft))]
dft.l=dft[1:(length(dft)-1)]

dst.t=dst[2:(length(dst))]
dst.l=dst[1:(length(dst)-1)]

vecm.1=lm(dft.t~ect.1)
summary(vecm.1)



vecm.2=lm(dst.t~ect.1)
summary(vecm.2)

#pre: when cointegration vector is know
#otherwise, we have to estimate it

vecm.1.pre=lm(dft.t~ect.pre)
vecm.2.pre=lm(dst.t~ect.pre)
summary(vecm.1.pre)
summary(vecm.2.pre)
####Newey-West Standard errors
library(sandwich)
sqrt(diag(vcovHC(vecm.1.pre)))
sqrt(diag(vcovHC(vecm.2.pre)))



# Example bivariate cointegrated system 
e1 <- rnorm(250,mean=0,sd=.5)
e2 <- rnorm(250,mean=0,sd=.5)
u.ar1 <- arima.sim(model=list(ar=.75),250,innov=e1)
y2 <- cumsum(e2)
y1 <- y2 + u.ar1
data1 <- cbind(y1,y2)

#ca.jo Johansen test
# trace statistic
test1 <- ca.jo(data1,ecdet="const",type="trace",K=2,spec="transitory")
  #ecdet
  #cointegration residual can be set to have one of the three properties:
    #Zero
    #constant
    #trend
  #need K to be atleast = 2 for function to work. The model requires a lag

##K Is the lag order in the VAR level
  #K refers to the lag length in the level of VAR model
summary(test1)

# max eigenvalue statistic
test2 <- ca.jo(data1,ecdet="const",type="eigen",K=2,spec="transitory")
  #instead of trace, we use eigen in type parameter.  
summary(test2)

#Johansen should be used to test how many cointegration vectors to specify
#not good for estimating model, it's very sensitive

# estimate restricted VECM (OLS regression of VECM)
model1 <- cajorls(test1,r=1)
summary(model1$rlm)
print(model1)





# Example trivariate cointegrated system (1 coint vector)
e1 <- rnorm(250,mean=0,sd=.5)
e2 <- rnorm(250,mean=0,sd=.5)
e3 <- rnorm(250,mean=0,sd=.5)
u.ar1 <- arima.sim(model=list(ar=.75),250,innov=e1)
y2 <- cumsum(e2)
y3 <- cumsum(e3)
y1 <- .5*y2 + .5*y3 + u.ar1
data2 <- cbind(y1,y2,y3)

# trace statistic
test1 <- ca.jo(data2,ecdet="const",type="trace",K=2,spec="transitory")
summary(test1)

# max eigenvalue statistic
test2 <- ca.jo(data2,ecdet="const",type="eigen",K=2,spec="transitory")
summary(test2)


# Example trivariate cointegrated system (2 coint vectors)
e1 <- rnorm(250,mean=0,sd=.5)
e2 <- rnorm(250,mean=0,sd=.5)
e3 <- rnorm(250,mean=0,sd=.5)
u.ar1 <- arima.sim(model=list(ar=.75),250,innov=e1)
v.ar1 <- arima.sim(model=list(ar=.75),250,innov=e2)
y3 <- cumsum(e3)
y1 <- y3 + u.ar1
y2 <- y3 + v.ar1
data3 <- cbind(y1,y2,y3)

# trace statistic
test1 <- ca.jo(data3,ecdet="const",type="trace",K=2,spec="transitory")
summary(test1)

# max eigenvalue statistic
test2 <- ca.jo(data3,ecdet="const",type="eigen",K=2,spec="transitory")
summary(test2)

# estimate restricted VECM
model1 <- cajorls(test1,r=2)
summary(model1$rlm)
print(model1)


# Example Zivot (2000)
# US-UK data 1976.02 - 1996.06
data <-read.table("usuk.txt",sep="",header=TRUE)

# get time series
st <- ts(data[,"USUKS"],start=c(1976,2),frequency=12)
ft <- ts(data[,"USUKF"],start=c(1976,2),frequency=12)
data <- cbind(st,ft)

# trace test
test1 <- ca.jo(data,ecdet="const",type="trace",K=2,spec="transitory")
summary(test1)

# max eigenvalue test
test2 <- ca.jo(data,ecdet="const",type="eigen",K=2,spec="transitory")
summary(test2)

# estimate restricted model
model1 <- cajorls(test1,r=1)
summary(model1$rlm)
print(model1)

