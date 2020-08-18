#########VAR AND UNIT ROOT TESTS######
##LIBRARIES vars AND urca####
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0327") # Set my working directory
library(vars)
library(urca)



# example Stock and Watson (2001)
data=read.table("SWdata.txt",sep="",header=TRUE)

# get time series
fedfr=data$ffrate[21:184]
inflr=400*(log(data$gdpd[21:184])-log(data$gdpd[20:183]))
unemr=data$urate[21:184]
gdpd=log(data$gdpd[21:184])#log-level of gdp deflator

# full sample
vardata0=ts(cbind(inflr,unemr,fedfr),start=1960,frequency=4)

# find optimal number of lags
info.crit=VARselect(vardata0,lag.max=8,type="const")
info.crit

# estimate VAR(2)
model0=VAR(vardata0,p=2,type="const")
summary(model0)

# get IRF
###Impact of shock to inflation
irf.inflr=irf(model0,impulse="inflr",n.ahead=24,ci=.9)
plot(irf.inflr)
###Impact of shock to unemployment

irf.unemr=irf(model0,impulse="unemr",n.ahead=24,ci=.9)
plot(irf.unemr)


irf.fedfr=irf(model0,impulse="fedfr",n.ahead=24,ci=.9)
plot(irf.fedfr)

# compute FEVD
fevd.var=fevd(model0,n.ahead=10)
plot(fevd.var)


#########SVAR IDENTIFICATION
amat=diag(3)

amat[2,1]=NA
amat[3,1]=NA
amat[3,2]=NA
amat
bmat=diag(3)
diag(bmat)=NA
bmat

svar0=SVAR(model0, estmethod = "direct", Amat = amat, Bmat = bmat,
hessian = TRUE, method = "BFGS")
summary(svar0)

############Overidentified Model####################
amat1=diag(3)

amat1
bmat1=diag(3)
diag(bmat1)=NA
bmat1
svar1=SVAR(model0, estmethod = "direct", Amat = amat1, Bmat = bmat1,
hessian = TRUE, method = "BFGS")
summary(svar1)
##svar is not identified if we have more structural parameters than reduced because there are more structural parameters than the
#reduced form parameters
##LR test compares the reduced form model with the restricted model
#LR=T(log(det(SigmaR))-det(sigmaU))sigmaR is the restricted variance-covariance matrix


###########Underidentified Model#############
amat2=diag(3)

amat2[1,2]=NA
amat2[3,1]=NA
amat2[3,2]=NA
amat2[2,1]=NA

bmat2=diag(3)
diag(bmat2)=NA
amat2
bmat2

svar2=SVAR(model0, estmethod = "direct", Amat = amat2, Bmat = bmat2,
           hessian = TRUE, method = "BFGS")
summary(svar2)


# get IRF
###Impact of shock to federal funds rate 
irf.fedfr.svar=irf(svar0,impulse="fedfr",n.ahead=24,ci=.9)
plot(irf.fedfr.svar)
plot(irf.fedfr)





######BQ Decomposition#############
dun=data$urate[21:184]-data$urate[20:183]

var.bq=ts(cbind(dun,fedfr),start=1960,frequency=4)
# find optimal number of lags
info.crit.bq=VARselect(var.bq,lag.max=10,type="const")
info.crit.bq

# estimate VAR(1)
model.bq=VAR(var.bq,p=1,type="const")
summary(model.bq)

svar.bq=BQ(model.bq)
summary(svar.bq)


# get IRF
irf.dun.bq=irf(svar.bq,impulse="dun",n.ahead=24,ci=.9,cumulative=FALSE)
irf.fedfr.bq=irf(svar.bq,impulse="fedfr",n.ahead=24,ci=.9,cumulative=FALSE)
irf.fedfr.bq

plot(irf.fedfr.bq)
plot(irf.dun.bq)



# get the impact of demand shocks on level of unemployment(unemr) 
# to get shocks to yt we need the cumulative sum of shocks to fedfr
# 
demand=cbind(cumsum(irf.fedfr.bq$irf$fedfr[,1]),
							irf.fedfr.bq$irf$fedfr[,2])

# plot IRFs - demand shocks
plot(demand[,1],type="l",col="black",lwd=2,ann=FALSE,
				xlim=c(0,24),ylim=c(-5,10))
lines(demand[,2],col="blue",lwd=2)
abline(h=0)
legend(x="topright",c("Unemployment response", "fed funds rate response"),
					col=c("black","blue"),lwd=2,cex=0.6)

