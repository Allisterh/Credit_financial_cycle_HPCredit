setwd("e:/835/prg_r") # Set my working directory
sink("assign1output.txt", append=TRUE, split=TRUE) #Save output from this script
lrgdp=read.csv("lrgdp.csv",header=T) # Load data
lrgdp=read.table("lrgdp.txt",header=T) # Load data
lrgdp.ts=ts(data=lrgdp,start=c(1947,1),frequency=4)
#Saving the plot
pdf(file = "graphsassign1.1.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(3,1))
## draw the plot
plot(lrgdp.ts,col="blue")
acf(lrgdp,main="acf of level of real GDP")
pacf(lrgdp,main="pacf of level of real GDP")


## close the device to do the drawing
dev.off()

#Estimate AR(1)
lrgdp.ar1=arima(lrgdp,order=c(1,0,0),method="ML")
lrgdp.ar1
#creating detrended log of real gdp. To do that, we need to run a regression of lrgdp
on a linear time trend

time=data.frame(c(1:length(lrgdp.ts))) #creating time dummy
names(time)[1] = "time"
lrgdp.df=data.frame(time,lrgdp)
level.mod=lm(lrgdp~time,data=lrgdp.df)
summary(level.mod)

dtlrgdp=resid(level.mod) #creating detrended log of real gdp which is the residual from the
dtlrgdp.ts=ts(data=dtlrgdp,start=c(1947,1),frequency=4)
#Saving the plot
pdf(file = "graphsassign1.2.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(3,1))
## draw the plot
plot(dtlrgdp.ts,col="blue")

acf(dtlrgdp,main="acf of detrended real GDP")

pacf(dtlrgdp,main="pacf of detrended real GDP")


## close the device to do the drawing
dev.off()


# Analysis with first difference of log of real GDP

dlrgdp=diff(lrgdp.ts)

#Saving the plot
pdf(file = "graphsassign1.3.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(3,1))
## draw the plot
plot(dlrgdp,col="blue")

acf(dlrgdp,main="acf of real GDP growth")

pacf(dlrgdp,main="pacf of real GDP growth")


## close the device to do the drawing
dev.off()

#Estimation of ARMA models for  dtlrgdp
model.dt.00=arima(dtlrgdp,order=c(0,0,0),method="ML")
model.dt.10=arima(dtlrgdp,order=c(1,0,0),method="ML")
model.dt.20=arima(dtlrgdp,order=c(2,0,0),method="ML")
model.dt.30=arima(dtlrgdp,order=c(3,0,0),method="ML",optim.control =
list(maxit = 2000))
model.dt.01=arima(dtlrgdp,order=c(0,0,1),method="ML")
model.dt.02=arima(dtlrgdp,order=c(0,0,2),method="ML")
model.dt.03=arima(dtlrgdp,order=c(0,0,3),method="ML")
model.dt.11=arima(dtlrgdp,order=c(1,0,1),method="ML")
model.dt.12=arima(dtlrgdp,order=c(1,0,2),method="ML")
model.dt.13=arima(dtlrgdp,order=c(1,0,3),method="ML")
model.dt.21=arima(dtlrgdp,order=c(2,0,1),method="ML")
model.dt.22=arima(dtlrgdp,order=c(2,0,2),method="ML")
model.dt.23=arima(dtlrgdp,order=c(2,0,3),method="ML")
model.dt.31=arima(dtlrgdp,order=c(3,0,1),method="ML")
model.dt.32=arima(dtlrgdp,order=c(3,0,2),method="ML",optim.control =
list(maxit = 2000))
model.dt.33=arima(dtlrgdp,order=c(3,0,3),method="ML",optim.control =
list(maxit = 2000))

#Estimation of ARMA models for dlrgdp
model.dl.00=arima(dlrgdp,order=c(0,0,0),method="ML")
model.dl.10=arima(dlrgdp,order=c(1,0,0),method="ML")
model.dl.20=arima(dlrgdp,order=c(2,0,0),method="ML")
model.dl.30=arima(dlrgdp,order=c(3,0,0),method="ML")
model.dl.01=arima(dlrgdp,order=c(0,0,1),method="ML")
model.dl.02=arima(dlrgdp,order=c(0,0,2),method="ML")
model.dl.03=arima(dlrgdp,order=c(0,0,3),method="ML")
model.dl.11=arima(dlrgdp,order=c(1,0,1),method="ML")
model.dl.12=arima(dlrgdp,order=c(1,0,2),method="ML")
model.dl.13=arima(dlrgdp,order=c(1,0,3),method="ML")
model.dl.21=arima(dlrgdp,order=c(2,0,1),method="ML")
model.dl.22=arima(dlrgdp,order=c(2,0,2),method="ML")
model.dl.23=arima(dlrgdp,order=c(2,0,3),method="ML")
model.dl.31=arima(dlrgdp,order=c(3,0,1),method="ML")
model.dl.32=arima(dlrgdp,order=c(3,0,2),method="ML")
model.dl.33=arima(dlrgdp,order=c(3,0,3),,optim.control =
list(maxit = 2000))


vec.aic.dt=vector()
vec.aic.dt[1]=AIC(model.dt.00)
vec.aic.dt[2]=AIC(model.dt.10)
vec.aic.dt[3]=AIC(model.dt.20)
vec.aic.dt[4]=AIC(model.dt.30)
vec.aic.dt[5]=AIC(model.dt.01)
vec.aic.dt[6]=AIC(model.dt.02)
vec.aic.dt[7]=AIC(model.dt.03)
vec.aic.dt[8]=AIC(model.dt.11)
vec.aic.dt[9]=AIC(model.dt.12)
vec.aic.dt[10]=AIC(model.dt.13)
vec.aic.dt[11]=AIC(model.dt.21)
vec.aic.dt[12]=AIC(model.dt.22)
vec.aic.dt[13]=AIC(model.dt.23)
vec.aic.dt[14]=AIC(model.dt.31)
vec.aic.dt[15]=AIC(model.dt.32)
vec.aic.dt[16]=AIC(model.dt.33)

vec.bic.dt=vector()
vec.bic.dt[1]=BIC(model.dt.00)
vec.bic.dt[2]=BIC(model.dt.10)
vec.bic.dt[3]=BIC(model.dt.20)
vec.bic.dt[4]=BIC(model.dt.30)
vec.bic.dt[5]=BIC(model.dt.01)
vec.bic.dt[6]=BIC(model.dt.02)
vec.bic.dt[7]=BIC(model.dt.03)
vec.bic.dt[8]=BIC(model.dt.11)
vec.bic.dt[9]=BIC(model.dt.12)
vec.bic.dt[10]=BIC(model.dt.13)
vec.bic.dt[11]=BIC(model.dt.21)
vec.bic.dt[12]=BIC(model.dt.22)
vec.bic.dt[13]=BIC(model.dt.23)
vec.bic.dt[14]=BIC(model.dt.31)
vec.bic.dt[15]=BIC(model.dt.32)
vec.bic.dt[16]=BIC(model.dt.33)


#AIC AND BIC FOR LOG DIFFERENCE MODEL

vec.aic.dl=vector()
vec.aic.dl[1]=AIC(model.dl.00)
vec.aic.dl[2]=AIC(model.dl.10)
vec.aic.dl[3]=AIC(model.dl.20)
vec.aic.dl[4]=AIC(model.dl.30)
vec.aic.dl[5]=AIC(model.dl.01)
vec.aic.dl[6]=AIC(model.dl.02)
vec.aic.dl[7]=AIC(model.dl.03)
vec.aic.dl[8]=AIC(model.dl.11)
vec.aic.dl[9]=AIC(model.dl.12)
vec.aic.dl[10]=AIC(model.dl.13)
vec.aic.dl[11]=AIC(model.dl.21)
vec.aic.dl[12]=AIC(model.dl.22)
vec.aic.dl[13]=AIC(model.dl.23)
vec.aic.dl[14]=AIC(model.dl.31)
vec.aic.dl[15]=AIC(model.dl.32)
vec.aic.dl[16]=AIC(model.dl.33)

vec.bic.dl=vector()
vec.bic.dl[1]=BIC(model.dl.00)
vec.bic.dl[2]=BIC(model.dl.10)
vec.bic.dl[3]=BIC(model.dl.20)
vec.bic.dl[4]=BIC(model.dl.30)
vec.bic.dl[5]=BIC(model.dl.01)
vec.bic.dl[6]=BIC(model.dl.02)
vec.bic.dl[7]=BIC(model.dl.03)
vec.bic.dl[8]=BIC(model.dl.11)
vec.bic.dl[9]=BIC(model.dl.12)
vec.bic.dl[10]=BIC(model.dl.13)
vec.bic.dl[11]=BIC(model.dl.21)
vec.bic.dl[12]=BIC(model.dl.22)
vec.bic.dl[13]=BIC(model.dl.23)
vec.bic.dl[14]=BIC(model.dl.31)
vec.bic.dl[15]=BIC(model.dl.32)
vec.bic.dl[16]=BIC(model.dl.33)


vec.aic.dt
vec.bic.dt


vec.aic.dl
vec.bic.dl


###Part (f) breaking the sample in two parts

dlrgdp.1=dlrgdp[1:139] ###1947:02-1981:04 (Note that one observation is lost in differencing
dlrgdp.2=dlrgdp[140:273]##1982:01-2015:02

dlrgdp.1.ts=ts(data=dlrgdp.1,start=c(1947,2),frequency=4)
dlrgdp.2.ts=ts(data=dlrgdp.2,start=c(1982,1),frequency=4)


######ACF AND PACF FOR THE FIRST SUB-SAMPLE
#Saving the plot
pdf(file = "graphsassign1.3.pdf")

par(mfrow = c(3,1))
## draw the plot
plot(dlrgdp.1.ts,col="blue")

acf(dlrgdp.1,main="acf of real GDP growth")

pacf(dlrgdp.1,main="pacf of real GDP growth")


## close the device to do the drawing
dev.off()

######ACF AND PACF FOR THE SECOND SUB-SAMPLE
#Saving the plot
pdf(file = "graphsassign1.4.pdf")

par(mfrow = c(3,1))
## draw the plot
plot(dlrgdp.2.ts,col="blue")

acf(dlrgdp.2,main="acf of real GDP growth")

pacf(dlrgdp.2,main="pacf of real GDP growth")


## close the device to do the drawing
dev.off()

auto.arima(dlrgdp.1.ts,ic="bic")
auto.arima(dlrgdp.2.ts,ic="bic")





