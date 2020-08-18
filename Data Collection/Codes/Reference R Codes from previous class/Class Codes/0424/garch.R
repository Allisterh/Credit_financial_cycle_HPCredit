setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0424") # Set my working directory

library(quantmod)
library(tseries)
library(forecast)

library(rugarch)
  #better assymetric garch than fgarch library


getSymbols("DEXUSUK",src="FRED") #Obtain Exchange rates from FRED

plot(DEXUSUK)
pp.test(na.omit(DEXUSUK))##Phillips-Perron Unit root test

DUSUK=diff(DEXUSUK)
DUSUK=na.omit(DUSUK)

pp.test(DUSUK)
plot(DUSUK)

Box.test(DUSUK, lag=12, type='Ljung')
##This suggest fitting an arma model to the 
##first difference of exchange rate
acf(DUSUK,lag=12)
pacf(DUSUK,lag=12)

#############ARCH/GARCH MODELS

model.ar1=arima(DUSUK,order=c(1,0,0),method="ML")
resid.ar1=resid(model.ar1)
Box.test(resid.ar1,lag=12,type='Ljung')

model.ar2=arima(DUSUK,order=c(2,0,0),method="ML")
resid.ar2=resid(model.ar2)
Box.test(resid.ar2,lag=12,type='Ljung')

Box.test(resid.ar1^2,lag=12,type='Ljung')
Box.test(resid.ar2^2,lag=12,type='Ljung')
  ##There is evidence of arch and garch affect
plot(resid.ar1^2)
  #pattern of clustering found here
options(digits=4)


arch10.spec = ugarchspec(variance.model = list(garchOrder=c(1,0)), 
                          mean.model = list(armaOrder=c(1,0)))
arch10.fit = ugarchfit(spec=arch10.spec, data=DUSUK, solver = "hybrid"
                             )
arch10.fit

    ##hybrid solver helps with converging problem, simulated MLE...


garch11.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                          mean.model = list(armaOrder=c(1,0)))
garch11.fit = ugarchfit(spec=garch11.spec, data=DUSUK,
                             )
garch11.fit
  ##for SME500 data, we see garch effect all the time

# Asymmetric garch
#

# Engle-Ng sign bias test
signbias(garch11.fit)
  ##insignificant, no sign bias here
  ##what does sign bias mean?

# Nelson's egarch model
egarch11.spec = ugarchspec(variance.model=list(model="eGARCH",
                                               garchOrder=c(1,1)),
                           mean.model=list(armaOrder=c(1,0)))
egarch11.fit = ugarchfit(egarch11.spec, DUSUK)
egarch11.fit
#the leverage effect parameter is alpha
  ##concur with the sign test bias

# GJR garch model(a version of TGARCH)
gjrgarch11.spec = ugarchspec(variance.model=list(model="gjrGARCH",
                                                 garchOrder=c(1,1)),
                             mean.model=list(armaOrder=c(1,0)))
gjrgarch11.fit = ugarchfit(gjrgarch11.spec, DUSUK)
gjrgarch11.fit
    ##gamma1 is leverage parameter

# aparch models
aparch11.1.spec = ugarchspec(variance.model=list(model="apARCH",
                                                 garchOrder=c(1,1)),
                             mean.model=list(armaOrder=c(1,0)),
                             fixed.pars=list(delta=1))

aparch11.1.fit = ugarchfit(aparch11.1.spec, DUSUK)
aparch11.1.fit
  ##gamma1 insignificant, what does it mean?
  
nic.garch11 = newsimpact(garch11.fit)
nic.egarch11 = newsimpact(egarch11.fit)
nic.gjrgarch11 = newsimpact(gjrgarch11.fit)
nic.aparch11.1 = newsimpact(aparch11.1.fit)

# compare information criteria
model.list = list(arch10=arch10.fit,
                  garch11 = garch11.fit,
                  egarch11 = egarch11.fit,
                  gjrgarch11 = gjrgarch11.fit,
                  aparch11.1 = aparch11.1.fit)
info.mat = sapply(model.list, infocriteria)
rownames(info.mat) = rownames(infocriteria(garch11.fit))
info.mat
  ##selection : garch11 has lowest values, arch10 and garch11 are close
  ##would still go with garch11 since it's the most popular model?

# show news impact curve from estimated garch(1,1) and egarch(1,1)
par(mfrow=c(2,2))
plot(nic.garch11$zx, type="l", lwd=2, col="blue", main="GARCH(1,1)", 
     nic.garch11$zy, ylab=nic.garch11$yexpr, xlab=nic.garch11$xexpr)
plot(nic.egarch11$zx, type="l", lwd=2, col="blue", main="EGARCH(1,1)", 
     nic.egarch11$zy, ylab=nic.egarch11$yexpr, xlab=nic.egarch11$xexpr)
plot(nic.gjrgarch11$zx, type="l", lwd=2, col="blue", main="TGARCH(1,1)", 
     nic.gjrgarch11$zy, ylab=nic.gjrgarch11$yexpr, xlab=nic.gjrgarch11$xexpr)
plot(nic.aparch11.1$zx, type="l", lwd=2, col="blue", main="APARCH(1,1,1)", 
     nic.aparch11.1$zy, ylab=nic.aparch11.1$yexpr, xlab=nic.aparch11.1$xexpr)

  ## plot of impact function
  ##which model to choose, use model selection criteria

# forecasts from competing models
#
garch11.fcst = ugarchforecast(garch11.fit, n.ahead=250)
egarch11.fcst = ugarchforecast(egarch11.fit, n.ahead=250)



