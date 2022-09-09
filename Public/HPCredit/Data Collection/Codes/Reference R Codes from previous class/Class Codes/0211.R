### Introduction ###
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes") # Set my working directory


###SIMULATING ARMA MODELS##########
# simulate gaussian white noise process
y=as.ts(rnorm(500))



# plot process
layout(1)
plot(y,col="blue",ylim=c(-3,3),ann=FALSE)
abline(h=0)

# Simulated AR(1) processes
y1=arima.sim(n=500,list(ar=.995,ma=0))#arima.sim is a function to simulate ARMA models
  #.995 is very close to stationary
y2=arima.sim(n=500,list(ar=.90,ma=0))
y3=arima.sim(n=500,list(ar=.75,ma=0))
y4=arima.sim(n=500,list(ar=.50,ma=0))
# plot time series
layout(matrix(c(1,2,3,4),2,2))
plot(y1,col="blue",ann=FALSE)
abline(h=0)
mtext("phi = 0.995")
plot(y2,col="blue",ann=FALSE)
abline(h=0)
mtext("phi = 0.90")
plot(y3,col="blue",ann=FALSE)
abline(h=0)
mtext("phi = 0.75")
plot(y4,col="blue",ann=FALSE)
abline(h=0)
mtext("phi= 0.50")
  #persistence shows on the graphs

acf(y1)
pacf(y1)
  #acf decline slowly, because .995^j decline slowly
acf(y4)
pacf(y4)
  #acf dropped much  faster since .5^j decreases faster

# Simulated AR(2) process
y.ar2=arima.sim(n=5000,list(ar=c(1.2,-0.4),ma=0),sd=0.5)

# plot time series
plot(y.ar2,col="blue")
abline(h=0)

# plot sample ACF and PACF 
acf(y.ar2)
pacf(y.ar2)

  #pacf shows negative value for second term

# estimate an ARMA(2,0)
model.ar2=arima(y.ar2,order=c(2,0,0),method="ML")
  #arima is the function to fit the data with, ML: maximum likelihood
model.ar2
  #this shows values for parameters of the estimated model
par(mar=c(1, 1, 1, 1))
par("mar")
tsdiag(model.ar2)
# find eigenvalues 
polyroot(c(-model.ar2$coef[2],-model.ar2$coef[1],1)) #will show imaginary values too
# test for autocorrelation in the residuals
Box.test(model.ar2$resid,lag=5,type="Ljung-Box")
Box.test(model.ar2$resid,lag=10,type="Ljung-Box")

###Estimate an AR(2) DGP by assuming that it's AR(1) and ARMA(1,1)
model.ar1=arima(y.ar2,order=c(1,0,0),method="ML")
model.arma11=arima(y.ar2,order=c(1,0,1),method="ML")

Box.test(model.ar1$resid,lag=5,type="Ljung-Box") #p-value is small -> wrong model
tsdiag(model.ar1)
# find eigenvalues
polyroot(c(-model.ar1$coef[1],1))
polyroot(c(1,-model.ar1$coef[1]))##roots of the characteristic equation


  #*** identify values of true model
vec.aic=vector() ###create an empty vector to store aic values

vec.aic[1]=AIC(model.ar2)
vec.aic[2]=AIC(model.ar1)
vec.aic[3]=AIC(model.arma11)

vec.aic #=> shows that aci for ar1 is the best

vec.bic=vector() ###create an empty vector to store bic values

vec.bic[1]=BIC(model.ar2)
vec.bic[2]=BIC(model.ar1)
vec.bic[3]=BIC(model.arma11)

vec.bic #=> again, this shows ar1 is the best

library(forecast)

#**auto.arima will find best model


auto.arima(y.ar2,ic="bic") #bic shows better/simler result, use it with a grain of salt
auto.arima(y.ar2) #without bic specification
####quantmod library for downloading data############
library(quantmod)   # Load the package


getSymbols("CPIAUCSL",src="FRED") #FRED good source for data
head(CPIAUCSL)
cpi.m.ts=ts(CPIAUCSL,start=c(1947,1),frequency=12)
plot(cpi.m.ts)
infl.m=(na.omit(diff(log(cpi.m.ts))))*1200 #take log and first difference
plot(infl.m)
cpi.q.ts=aggregate(cpi.m.ts,nfrequency=4,mean) #aggregate monthly to quarterly data

infl.q.ts=(na.omit(diff(log(cpi.q.ts))))*400
plot(infl.q.ts)
acf(infl.q.ts)
acf(coredata(infl.q.ts)) #*get around the problem of data visualization
pacf(coredata(infl.q.ts))
auto.arima(infl.q.ts,ic="bic")


  #next is unemployment data
par (mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
getSymbols("UNRATE",src="FRED",from="1960-01", to="2018-12")  # Download data from FRED.
u=UNRATE
plot(u)
acf(u)
pacf(u)
u.m.ts=ts(u,start=c(1948,1),frequency=12)
u.q.ts=aggregate(u.m.ts,nfrequency=4,mean)
acf(coredata(u.m.ts))
pacf(coredata(u.m.ts))
du.m=na.omit(diff(u.m.ts))
du.q=na.omit(diff(u.q.ts))
plot(du.q)

acf(coredata(du.m))
pacf(coredata(du.m))
acf(du.m)
acf(coredata(du.q))
pacf(coredata(du.q))




