#########VAR MODELS######
##LIBRARIES vars AND urca####
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0313/R_Scripts_for_VAR_Models") # Set my working directory
library(vars)




# example Stock and Watson (2001)
data=read.table("SWdata.txt",sep="",header=TRUE)

# get time series
fedfr=data$ffrate[21:184]
inflr=400*(log(data$gdpd[21:184])-log(data$gdpd[20:183]))
unemr=data$urate[21:184]
gdpd=log(data$gdpd[21:184])#log-level of gdp deflator

# full sample
vardata0=ts(cbind(inflr,unemr,fedfr),start=1960,frequency=4)
# plot processes
layout(1)
plot(vardata0,type="l",col="blue",main="Stock and Watson (2001) Data")



# find optimal number of lags
info.crit=VARselect(vardata0,lag.max=8,type="const")
info.crit

# estimate VAR(2)
model0=VAR(vardata0,p=2,type="const")
model1=VAR(vardata0,ic="SC",type="const")
summary(model0)
summary(model1)
# get roots (modulus)
roots(model0)
  #since all values <1 -> this means the series is stationary
# model evaluation
var2.serial=serial.test(model0,lags.pt=5,type="PT.adjusted")
var2.serial
  # test for serial correlation
  #p values cloes to 0 is not a good news. we left something out here

# test for Granger causality
causality(model0,cause=c("inflr","unemr"))$Granger
  #null hypo H0: inflr unemr do not Granger-cause fedfr
causality(model0,cause=c("inflr"))$Granger
  #without second parameter, the function gives causality for all other variables
  #p=0 means lags from other variables are important too
causality(model0,cause=c("fedfr","unemr"))$Granger
causality(model0,cause=c("inflr","fedfr"))$Granger



# get IRF (impulse responses)
###Impact of shock to inflation
irf.inflr=irf(model0,impulse="inflr",n.ahead=24,ci=.9) #impact of shocks to inflation
plot(irf.inflr)
  #apply cholesky decomposition, lower triangular matrix
  #when inflation goes up, fed fundrate (i) goes up accordingly and this would
    #induce a contractionary effect -> unemployment goes up

###Impact of shock to unemployment

irf.unemr=irf(model0,impulse="unemr",n.ahead=24,ci=.9)
plot(irf.unemr)

###Impact of shock to federal funds rate
  #unexpected increase in fed fund rate
irf.fedfr=irf(model0,impulse="fedfr",n.ahead=24,ci=.9)
plot(irf.fedfr)

# compute FEVD (forecast error decomposition) useful graphs showing contribution
  #contribution to variances in 3 variables
fevd.var=fevd(model0,n.ahead=16)
plot(fevd.var)


# set up data set **now ordering matters**
vardata1=ts(cbind(unemr,inflr,fedfr),start=1960,frequency=4)

# estimate VAR(2)  unemr->inflr->fedfr
model1=VAR(vardata1,p=2,type="const")
summary(model1)

# get IRF
irf.fedfr=irf(model1,impulse="fedfr",n.ahead=24,ci=.9)
plot(irf.fedfr)

irf.inflr=irf(model1,impulse="inflr",n.ahead=24,ci=.9)
plot(irf.inflr)

irf.unemr=irf(model1,impulse="unemr",n.ahead=24,ci=.9)
plot(irf.unemr)


# set up data set **now ordering matters**
vardata2=ts(cbind(fedfr,unemr,inflr),start=1960,frequency=4)

# estimate VAR(2)  fedfr->unemr->inflr
model2=VAR(vardata2,p=2,type="const")
summary(model2)

# get IRF
irf.fedfr.2=irf(model2,impulse="fedfr",n.ahead=24,ci=.9)
plot(irf.fedfr.2)

irf.inflr.2=irf(model2,impulse="inflr",n.ahead=24,ci=.9)
plot(irf.inflr.2)

irf.unemr.2=irf(model2,impulse="unemr",n.ahead=24,ci=.9)
plot(irf.unemr.2)


