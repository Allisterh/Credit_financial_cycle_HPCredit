## 1. Load data ----
rm(list=ls())
library(vars)
library(stats)

country = 'GB'
#setwd("D:/GitHub/HPCredit/Data Collection/1.Latest")

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#-------------
## Data input
### Data source values available from 1988-07-01 -> 2019-10-01

filepath = sprintf("../../../Data Collection/1.Latest/MergedData_%s.txt",country)
df <- read.table(filepath, header=TRUE, sep=",")

# Limit series data to after 1990 ---------


#Cycles var name list
varlist = c("ID", "date", "Credit", "HPIndex")
df1 = df[varlist]

df1$Credit <- 100*log(df1$Credit)
df1$HPIndex <- 100*log(df1$HPIndex)


## filter Credit and HPIndex series using HP filter @lambda =400000
df1$Credit_HPtrend = mFilter::hpfilter(df1$Credit, type = "lambda", freq = 26000)$trend
df1$Credit_HPcycle = mFilter::hpfilter(df1$Credit, type = "lambda", freq = 26000)$cycle
df1$HPIndex_HPtrend = mFilter::hpfilter(df1$HPIndex, type = "lambda", freq = 26000)$trend
df1$HPIndex_HPcycle = mFilter::hpfilter(df1$HPIndex, type = "lambda", freq = 26000)$cycle


df2<-cbind(df['Credit'], df['HPIndex'], df1['Credit_HPcycle'],
           df1['HPIndex_HPcycle'], df1['Credit_HPtrend'], df1['HPIndex_HPtrend'])


filepath = sprintf("../../../Data Collection/1.Latest/MergedData_Matlab_%s.txt",country)
write.table(df2,filepath, sep=",")


var.aic <- VAR(df1[c('Credit_HPcycle','HPIndex_HPcycle')], type = "none", lag.max = 2, ic = "AIC") #requires different package, refer to class codes

#var.aic <- VAR(df, type = "none")
summary(var.aic)
var.aic$varresult$Credit_HPcycle$coefficients
var.aic$varresult$HPIndex_HPcycle$coefficients

p1<-summary(var.aic)[["varresult"]][["Credit_HPcycle"]][["coefficients"]][,1]
p1<-as.numeric(p1)
p11<-p1
p11[2] <- p1[3]
p11[3] <- p1[2]
p1<-p11
p2<-summary(var.aic)[["varresult"]][["HPIndex_HPcycle"]][["coefficients"]][,1]
p2<-as.numeric(p2)
p22<-p2
p22[2] <- p2[3]
p22[3] <- p2[2]
p2<-p22
p3<-summary(var.aic)[["varresult"]][["Credit_HPcycle"]][["sigma"]]^2
p4<-summary(var.aic)[["varresult"]][["HPIndex_HPcycle"]][["sigma"]]^2

prior<- c(p1,p2,p3,p4)

prior_path = sprintf("../Priors/prior_VAR2x_%s.txt",country)
write.csv(prior, prior_path)

# from uc_yc.R v2 code

############ Part 2: AR2
#AR2 Regression code
#Export Priors from separate AR(2) of HP-filter

y = df1[,"Credit_HPcycle"]
h = df1[,"HPIndex_HPcycle"]
cov1 = cov(y,h)
# araic = ar.mle(x, aic = TRUE, order.max = NULL, na.action = na.fail,
#        demean = TRUE)
# summary(araic)

library(forecast)
library(stats)
auto.arima(y)
auto.arima(h)

arima1 = arima(y,order=c(2,0,0),method="ML", include.mean = FALSE)
c1=arima1$coef
#c2=sqrt(diag(vcov(arima1)))
c5=arima1$sigma2
arima1 = arima(h,order=c(2,0,0),method="ML", include.mean = FALSE)
c3=arima1$coef
#c4=sqrt(diag(vcov(arima1)))
c6=arima1$sigma2

prior = cbind(t(c1),t(c3),c5,c6)
prior = as.numeric(prior)
  #prior structure: 8 in total
  #y: ar1 ar2, h:ar1 ar2
  #4 respective standard error
  #2 sigma^2 (residual sum of square)
#combine 
#export to prior_VAR2_US.txt

prior_path = sprintf("../Priors/prior_VAR2_%s.txt",country)
write.csv(prior, prior_path)


# library(PerformanceAnalytics)
# StdDev(y)
# StdDev(h)
# sd(y)


#### Part 3 Std Dev of trends ==============
#Calculate StdDev of trends

#df1 = subset(df, date > as.Date("1988-12-31"))

#Cycles var name list
#varlist2 = c("ID", "date", "Credit_HPtrend", "HPIndex_HPtrend")
#df1 = df[varlist2]
#var.aic <- VAR(df1[c(3,4)], type = "none", lag.max = 5, ic = "AIC") #requires different package, refer to class codes
#summary(var.aic)

y = df1[,"Credit_HPtrend"]
h = df1[,"HPIndex_HPtrend"]

cov2 = cov(y,h)
cov_yh = c(cov1,cov2)
prior_path = sprintf("../Priors/prior_corr_%s.txt",country)
write.csv(cov_yh, prior_path)

library(forecast)
auto.arima(y)
auto.arima(h)

arima1 = arima(y,order=c(0,1,0),method="ML", include.mean = FALSE)
c1=arima1$coef
c2=sqrt(arima1$sigma2)

arima1 = arima(h,order=c(0,1,0),method="ML", include.mean = FALSE)
c3=arima1$coef
c4=sqrt(arima1$sigma2)

prior = cbind(c2,c4)
prior = as.numeric(prior)

prior_path = sprintf("../Priors/prior_trend_%s.txt",country)
write.csv(prior, prior_path)

# x = seq(1, length(y),1)
# f <- lm(y~x)
# summary(f) #estimate drift coefficient
# cor(y,h)

# araic = ar.mle(x, aic = TRUE, order.max = NULL, na.action = na.fail,
#        demean = TRUE)
# summary(araic)

library(forecast)
auto.arima(y)
auto.arima(h)
# 
# StdDev(y)
# StdDev(h)
# sd


