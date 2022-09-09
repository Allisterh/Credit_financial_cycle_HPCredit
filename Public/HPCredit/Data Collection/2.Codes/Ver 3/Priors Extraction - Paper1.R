## 1. Load data ----
rm(list = ls())
library(vars)
library(stats)
library(forecast)


library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
setwd("../../1.Latest/Paper1")

source("../../2.Codes/Ver 3/HPfilters/OneSidedHPfilterfunc.R")
# source("../../2.Codes/Ver 3/HPfilters/OneSidedSTMfilterfunc.R")

#-------------
filepath <- "shortlistofCountries.csv"
name1 <- read.table(filepath, header = TRUE, sep = ",")

for (i in seq_len(nrow(name1))) {
country <- name1[i, ]


filepath <- sprintf("creditHPI_%s.csv", country)
df <- read.table(filepath, header = TRUE, sep = ",")

# Log transformation
df$creditlog <- ts(100 * log(df$Credit))
df$hpilog <- ts(100 * log(df$HPI))

df$credit_trend <- filterHP(df$creditlog, lambda = 125000)[, "trend"]
df$credit_cycle <- df$creditlog - df$credit_trend

df$hpi_trend <- filterHP(df$hpilog, lambda = 125000)[, "trend"]
df$hpi_cycle <- df$hpilog - df$hpi_trend


filepath <- sprintf("creditHPI_Matlab_%s.csv", country)
write.table(df[,-c(1:2)], filepath, sep = ",", row.names = FALSE, col.names=FALSE)


#Cycles var name list
varlist <- c("ID", "date", "credit_cycle", "hpi_cycle")
df1 <- df[varlist]

var_aic <- VAR(df1[c(3,4)], type = "none", lag.max = 2, ic = "AIC")
#requires different package, refer to class codes
#var_aic <- VAR(df, type = "none")
summary(var_aic)
var_aic$varresult$credit_cycle$coefficients
var_aic$varresult$hpi_cycle$coefficients

summary(var_aic)$covres
summary(var_aic)$corres

p1 <- summary(var_aic)[["varresult"]][["credit_cycle"]][["coefficients"]][, 1]
p1 <- as.numeric(p1)
p11 <- p1
p11[2] <- p1[3]
p11[3] <- p1[2]
p1 <- p11
p2 <- summary(var_aic)[["varresult"]][["hpi_cycle"]][["coefficients"]][, 1]
p2 <- as.numeric(p2)
p22 <- p2
p22[2] <- p2[3]
p22[3] <- p2[2]
p2 <- p22
p3 <- summary(var_aic)[["varresult"]][["credit_cycle"]][["coefficients"]][, 2]
p3 <- as.numeric(p3)
p33 <- p3
p33[2] <- p3[3]
p33[3] <- p3[2]
p3 <- p33
p4 <- summary(var_aic)[["varresult"]][["hpi_cycle"]][["coefficients"]][, 2]
p4 <- as.numeric(p4)
p44 <- p4
p44[2] <- p4[3]
p44[3] <- p4[2]
p4 <- p44
p5 <- summary(var_aic)[["varresult"]][["credit_cycle"]][["sigma"]]
p6 <- summary(var_aic)[["varresult"]][["hpi_cycle"]][["sigma"]]
p7 <- summary(var_aic)[["corres"]][[2]]
# https://stackoverflow.com/questions/36690093/extract-coefficients-and-variance-covariance-matrix-from-var-output-estimated-w

prior <- c(p1, p2, p3, p4, p5, p6, p7)
prior_path <- sprintf("Priors/prior_VAR2x_%s.csv", country)
write.csv(prior, prior_path)

# from uc_yc.R v2 code


############ Part 2: AR2
#AR2 Regression code
#Export Priors from separate AR(2) of HP-filter

y = df1[, c(3)]
h = df1[, c(4)]
cor1 = cov(y, h)
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
c3=sqrt(arima1$sigma2)
arima1 = arima(h,order=c(2,0,0),method="ML", include.mean = FALSE)
c4=arima1$coef
c5=sqrt(diag(vcov(arima1)))
c6=sqrt(arima1$sigma2)
prior = cbind(c1,c4,c2,c5) 
prior = as.numeric(prior)
prior = c(prior, c3, c6)
  #prior structure: 8 in total
  #y: ar1 ar2, h:ar1 ar2
  #4 respective standard error 
#combine 
#export to prior_VAR2_US.txt

prior_path = sprintf("Priors/prior_VAR2_%s.csv",country)
write.csv(prior, prior_path)

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

#df1 = subset(df, date > as.Date("1988-12-31"))

#Cycles var name list
varlist2 = c("ID", "date", "credit_trend", "hpi_trend")
df1 = df[varlist2]
var_aic <- VAR(df1[c(3,4)], type = "both", lag.max = 1, ic = "AIC") #requires different package, refer to class codes
summary(var_aic)

y = df1[,c(3)]
h = df1[,c(4)]

# cor2 = cor(y,h)
# cor_yh = c(cor1,cor2)
# prior_path = sprintf("../1.2.Priors/prior_corr_%s.txt",country)
# write.csv(cor_yh, prior_path)

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

prior_path = sprintf("Priors/prior_trend_%s.csv",country)
write.csv(prior, prior_path)

# x = seq(1, length(y),1)
# f <- lm(y~x)
# summary(f) #estimate drift coefficient
# cor(y,h)

# araic = ar.mle(x, aic = TRUE, order.max = NULL, na.action = na.fail,
#        demean = TRUE)
# summary(araic)

# library(forecast)
# auto.arima(y)
# auto.arima(h)
# 
# StdDev(y)
# StdDev(h)
# sd
}