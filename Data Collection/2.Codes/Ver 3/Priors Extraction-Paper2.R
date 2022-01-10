## 1. Load data ----
rm(list=ls())
library(vars)
library(stats)

country = 'JP'

# Start date is 1989-01-01
startdate_uc = '1988-07-01' #push back 2 periods to extract priors for VAR(2) process
startdate_var = '1989-01-01'
# End date is pre-covid 2020-01-01
enddate = '2021-04-01'

setwd(dirname(getActiveDocumentContext()$path))
setwd("../../1.Latest/Paper2")


#-------------

filepath = sprintf("MergedData_%s.txt",country)
df <- read.table(filepath, header=TRUE, sep=",")

# HP filter
df$HPIndex_HPtrend = mFilter::hpfilter(df$HPIndex, type = "lambda", freq = 1600)$trend
df$HPIndex_HPcycle = mFilter::hpfilter(df$HPIndex, type = "lambda", freq = 1600)$cycle
df$Credit_HPtrend = mFilter::hpfilter(df$credit, type = "lambda", freq = 1600)$trend
df$Credit_HPcycle = mFilter::hpfilter(df$credit, type = "lambda", freq = 1600)$cycle


df$date=as.Date(df$date)
df_matlab = subset(df, date >= as.Date(startdate_uc))
df_matlab = subset(df, date <= as.Date(enddate))
varlist = c("credit", "HPIndex", "Credit_HPcycle", "HPIndex_HPcycle", "Credit_HPtrend", "HPIndex_HPtrend")
df_matlab <-df_matlab[varlist]
filepath = sprintf("MergedData_Matlab_%s.txt",country)
write.table(df_matlab, filepath, sep=",")

# Limit series data to after 1990 ---------
df1 = subset(df, date >= as.Date(startdate_var))

#Cycles var name list
varlist = c("ID", "date", "Credit_HPcycle", "HPIndex_HPcycle")
df1 = df1[varlist]

var.aic <- VAR(df1[c("Credit_HPcycle","HPIndex_HPcycle")], type = "none", lag.max = 2, ic = "AIC") #requires different package, refer to class codes
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
p3<-summary(var.aic)[["varresult"]][["Credit_HPcycle"]][["sigma"]]
p4<-summary(var.aic)[["varresult"]][["HPIndex_HPcycle"]][["sigma"]]

prior<- c(p1,p2,p3,p4)

prior_path = sprintf("../../1.2.Priors/prior_VAR2x_%s.txt",country)
write.csv(prior, prior_path)

# from uc_yc.R v2 code

############ Part 2: AR2
#AR2 Regression code
#Export Priors from separate AR(2) of HP-filter

y = df1[,"Credit_HPcycle"]
h = df1[,"HPIndex_HPcycle"]
cor1 = cov(y,h)
# araic = ar.mle(x, aic = TRUE, order.max = NULL, na.action = na.fail,
#        demean = TRUE)
# summary(araic)

library(forecast)
library(stats)
auto.arima(y)
auto.arima(h)

arima1 = arima(y,order=c(2,0,0),method="ML", include.mean = FALSE)
c1=arima1$coef
# c2=sqrt(diag(vcov(arima1)))
c2=arima1$sigma2
arima1 = arima(h,order=c(2,0,0),method="ML", include.mean = FALSE)
c3=arima1$coef
# c4=sqrt(diag(vcov(arima1)))
c4=arima1$sigma2

prior = cbind(c1,c3,c2,c4)
prior = as.numeric(prior)
  #prior structure: 8 in total
  #y: ar1 ar2, h:ar1 ar2
  #4 respective standard error 
#combine 
#export to prior_VAR2_US.txt

prior_path = sprintf("../../1.2.Priors/prior_VAR2_%s.txt",country)
write.csv(prior, prior_path)


prior = cbind(c1,c2)
prior = as.numeric(prior)
prior_path = sprintf("../../1.2.Priors/prior_VAR2_credit_%s.txt",country)
write.csv(prior, prior_path)

prior = cbind(c3,c4)
prior = as.numeric(prior)
prior_path = sprintf("../../1.2.Priors/prior_VAR2_hpi_%s.txt",country)
write.csv(prior, prior_path)


# library(PerformanceAnalytics)
# StdDev(y)
# StdDev(h)
# sd(y)


#### Part 3 Std Dev of trends ==============
#Calculate StdDev of trends

df1 = subset(df, date >= as.Date(startdate_var))

#Cycles var name list
varlist2 = c("ID", "date", "Credit_HPtrend", "HPIndex_HPtrend")
df1 = df[varlist2]
var.aic <- VAR(df1[c(3,4)], type = "none", lag.max = 5, ic = "AIC") #requires different package, refer to class codes
summary(var.aic)

y = df1[,c(3)]
h = df1[,c(4)]

cor2 = cor(y,h)
cor_yh = c(cor1,cor2)
prior_path = sprintf("../../1.2.Priors/prior_corr_%s.txt",country)
write.csv(cor_yh, prior_path)

library(forecast)
auto.arima(y)
auto.arima(h)

arima1 = arima(y,order=c(1,0,0),method="ML", include.mean = FALSE)
c1=arima1$coef
c2=sqrt(arima1$sigma2)

arima1 = arima(h,order=c(1,0,0),method="ML", include.mean = FALSE)
c3=arima1$coef
c4=sqrt(arima1$sigma2)

prior = cbind(c2,c4)
prior = as.numeric(prior)

prior_path = sprintf("../../1.2.Priors/prior_trend_%s.txt",country)
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


