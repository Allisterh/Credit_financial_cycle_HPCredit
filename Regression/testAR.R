# set.seed(123) # Reset random number generator for reasons of reproducability
# 
# # Generate sample
# t <- 200 # Number of time series observations
# k <- 2 # Number of endogenous variables
# p <- 2 # Number of lags
# 
# # Generate coefficient matrices
# A.1 <- matrix(c(-.3, .6, -.4, .5), k) # Coefficient matrix of lag 1
# A.2 <- matrix(c(-.1, -.2, .1, .05), k) # Coefficient matrix of lag 2
# A <- cbind(A.1, A.2) # Companion form of the coefficient matrices
# 
# # Generate series
# series <- matrix(0, k, t + 2*p) # Raw series with zeros
# for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,0.5)
#   series[, i] <- A.1%*%series[, i-1] + A.2%*%series[, i-2] + rnorm(k, 0, .5)
# }
# 
# series <- ts(t(series[, -(1:p)])) # Convert to time series format
# names <- c("V1", "V2") # Rename variables
# 
# plot.ts(series) # Plot the series
# 
# library(vars) # Load package
# 
# var.1 <- VAR(series, 2, type = "none") # Estimate the model
# 
# var.aic <- VAR(series, type = "none", lag.max = 5, ic = "AIC")
# summary(var.aic)
########## Part 1: VAR2x
rm(list=ls())
library(vars)
library(stats)

# VAR estimation ========

country = 'US'
working_dir=("D:/GitHub/HPCredit/Regression/")
setwd(working_dir) 


#-------------

HP_filepath = sprintf("../Data/Input/HPindex_HPfilter_%s.txt",country)
df1 <- read.table(HP_filepath, header=TRUE, sep=",")

Credit_filepath = sprintf("../Data/Input/Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy

# Limit series data to after 1990 ---------
df <- merge(df1, df2, by=c("ID","date"))
df <- na.omit(df)
df = subset(df, date > as.Date("1989-12-31"))
varlist = c("Credit", "HPIndex")
df3 = df[varlist]
filepath = sprintf("../Data/Input/data_trimmed_%s.txt",country)
write.csv(df3, filepath)

df <- merge(df1, df2, by=c("ID","date"))
df <- na.omit(df)
#Cycles var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "HPIndex_HPcycle")
df = df[varlist2]
df = subset(df, date > as.Date("1989-09-30")) #limit data to after 1990
varlist2 = c("Credit_HPcycle","HPIndex_HPcycle")
df3 = df[varlist2]

file_path = sprintf("../Data/Input/Cycles_trimmed_%s.txt",country)
write.csv(df3, file_path)

var.aic <- VAR(df[c(3,4)], type = "none", lag.max = 2, ic = "AIC") #requires different package, refer to class codes
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

prior_path = sprintf("Priors/prior_VAR2x_%s.txt",country)
write.csv(prior, prior_path)

# from uc_yc.R v2 code

############ Part 2: AR2
#AR2 Regression code
#Export Priors from separate AR(2) of HP-filter

y = df[,c(3)]
h = df[,c(4)]
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
c2=sqrt(diag(vcov(arima1)))
arima1 = arima(h,order=c(2,0,0),method="ML", include.mean = FALSE)
c3=arima1$coef
c4=sqrt(diag(vcov(arima1)))

prior = cbind(c1,c3,c2,c4)
prior = as.numeric(prior)
  #prior structure: 8 in total
  #y: ar1 ar2, h:ar1 ar2
  #4 respective standard error 
#combine 
#export to prior_VAR2_US.txt

prior_path = sprintf("Priors/prior_VAR2_%s.txt",country)
write.csv(prior, prior_path)


# library(PerformanceAnalytics)
# StdDev(y)
# StdDev(h)
# sd(y)


#### Part 3 Std Dev of trends ==============
#Calculate StdDev of trends

df <- merge(df1, df2, by=c("ID","date"))
df <- na.omit(df)
df = subset(df, date > as.Date("1989-12-31"))

#Cycles var name list
varlist2 = c("ID", "date", "Credit_HPtrend", "HPIndex_HPtrend")
df = df[varlist2]
var.aic <- VAR(df[c(3,4)], type = "none", lag.max = 5, ic = "AIC") #requires different package, refer to class codes
summary(var.aic)

y = df[,c(3)]
h = df[,c(4)]

cor2 = cor(y,h)
cor_yh = c(cor1,cor2)
prior_path = sprintf("Priors/prior_corr_%s.txt",country)
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

prior_path = sprintf("Priors/prior_trend_%s.txt",country)
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


