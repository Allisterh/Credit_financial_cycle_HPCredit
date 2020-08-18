set.seed(123) # Reset random number generator for reasons of reproducability

# Generate sample
t <- 200 # Number of time series observations
k <- 2 # Number of endogenous variables
p <- 2 # Number of lags

# Generate coefficient matrices
A.1 <- matrix(c(-.3, .6, -.4, .5), k) # Coefficient matrix of lag 1
A.2 <- matrix(c(-.1, -.2, .1, .05), k) # Coefficient matrix of lag 2
A <- cbind(A.1, A.2) # Companion form of the coefficient matrices

# Generate series
series <- matrix(0, k, t + 2*p) # Raw series with zeros
for (i in (p + 1):(t + 2*p)){ # Generate series with e ~ N(0,0.5)
  series[, i] <- A.1%*%series[, i-1] + A.2%*%series[, i-2] + rnorm(k, 0, .5)
}

series <- ts(t(series[, -(1:p)])) # Convert to time series format
names <- c("V1", "V2") # Rename variables

plot.ts(series) # Plot the series

library(vars) # Load package

var.1 <- VAR(series, 2, type = "none") # Estimate the model

var.aic <- VAR(series, type = "none", lag.max = 5, ic = "AIC")
summary(var.aic)

# VAR estimation ========
country = 'US'
working_dir=sprintf("D:/GitHub/HPCredit/Data/Input/HPindex_HPfilter_%s.txt",country)

Credit_filepath = sprintf("../../../Data/Input/Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy

df <- merge(df1, df2, by=c("ID","date"))
df <- na.omit(df)
#Cycles var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "HPIndex_HPcycle")
df = df[varlist2]

var.aic <- VAR(df, type = "none", lag.max = 5, ic = "AIC")
summary(var.aic)

x = df[,c(3)]
y = df[,c(4)]

araic = ar.mle(x, aic = TRUE, order.max = NULL, na.action = na.fail,
       demean = TRUE)
summary(araic)

library(forecast)
auto.arima(x)
arima(y,order=c(2,0,0),method="ML")