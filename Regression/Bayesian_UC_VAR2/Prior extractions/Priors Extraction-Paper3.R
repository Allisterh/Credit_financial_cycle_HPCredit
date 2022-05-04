## 1. Load data ----
rm(list=ls())
library(vars)
library(stats)
library(mFilter)


setwd(dirname(getActiveDocumentContext()$path))
setwd("../../1.Latest/Paper3")

#country = 'AR'

filepath = "../../../../3rdPaper/Data/Processed/countrylist.txt"
countrylist <- read.table(filepath, header=TRUE, sep=",")




# Start date is 1989-01-01
startdate_uc = '1985-07-01' #push back 2 periods to extract priors for VAR(2) process
startdate_var = '1985-01-01'
# End date is pre-covid 2020-01-01
enddate = '2018-10-01'




#-------------

filepath = "credit.csv"
df0 <- read.table(filepath, header=TRUE, sep=",")

for(i in 1:nrow(countrylist)){

country = countrylist[i,]
  
df <- df0[c('date',country)]
names(df) <- c('date','credit')

# Log transformation
#df$creditlog <- 100*log(df$credit)

# HP filter
df <- na.omit(df)

df$Credit_HPtrend = mFilter::hpfilter(df$credit, type = "lambda", freq = 400000)$trend
df$Credit_HPcycle = mFilter::hpfilter(df$credit, type = "lambda", freq = 400000)$cycle


df$date=as.Date(df$date)
df_matlab = subset(df, date >= as.Date(startdate_var))
df_matlab = subset(df_matlab, date <= as.Date(enddate))
varlist = c("credit", "Credit_HPcycle", "Credit_HPtrend")
df_matlab <-df_matlab[varlist]
filepath = sprintf("MergedData_Matlab_%s.txt",country)
write.table(df_matlab, filepath, sep=",")

# Limit series data to after 1990 ---------
df1 = subset(df, date >= as.Date(startdate_var) & date <= as.Date(enddate))

#Cycles var name list
varlist = c("date", "Credit_HPcycle")
df1 = df1[varlist]

# from uc_yc.R v2 code

############ Part 2: AR2
#AR2 Regression code
#Export Priors from separate AR(2) of HP-filter

y = df1[,"Credit_HPcycle"]
# araic = ar.mle(x, aic = TRUE, order.max = NULL, na.action = na.fail,
#        demean = TRUE)
# summary(araic)

library(forecast)
library(stats)
#auto.arima(y)

arima1 = arima(y,order=c(2,0,0),method="ML", include.mean = FALSE)
c1=arima1$coef

# c2=sqrt(diag(vcov(arima1)))
c2=sqrt(arima1$sigma2)


  #prior structure: 8 in total
  #y: ar1 ar2, h:ar1 ar2
  #4 respective standard error 
#combine 
#export to prior_VAR2_US.txt
c1 = t(c1)
prior = cbind(c1,c2)
prior = as.numeric(prior)
prior_path = sprintf("../../1.Latest/Paper3/Priors/prior_VAR2_credit_%s.txt",country)
write.csv(prior, prior_path)

# library(PerformanceAnalytics)
# StdDev(y)
# StdDev(h)
# sd(y)


#### Part 3 Std Dev of trends ==============
#Calculate StdDev of trends

#Cycles var name list
varlist2 = c("date", "Credit_HPtrend")
df1 = df[varlist2]

y = df1["Credit_HPtrend"]

library(forecast)
auto.arima(y)

arima1 = arima(y,order=c(0,0,0),method="ML", include.mean = FALSE)
c1=arima1$coef
c2=sqrt(arima1$sigma2)

prior = c2
prior = as.numeric(prior)

prior_path = sprintf("../../1.Latest/Paper3/Priors/prior_trend_%s.txt",country)
write.csv(prior, prior_path)
}
# x = seq(1, length(y),1)
# f <- lm(y~x)
# summary(f) #estimate drift coefficient
# cor(y,h)

# araic = ar.mle(x, aic = TRUE, order.max = NULL, na.action = na.fail,
#        demean = TRUE)
# summary(araic)

#library(forecast)
#auto.arima(y)
# 
# StdDev(y)
# StdDev(h)
# sd


