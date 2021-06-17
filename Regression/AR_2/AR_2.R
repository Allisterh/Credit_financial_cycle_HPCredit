library(zoo)
library(rucm)


rm(list=ls())


country = 'US'
setwd("D:/GitHub/HPCredit/Data Collection/1.Latest")

#-------------

filepath = sprintf("MergedData_%s.txt",country)
df <- read.table(filepath, header=TRUE, sep=",")


# # Limit series data to after 1990 ---------
# 
#Cycles var name list
varlist = c("date", "Credit")
df1 = df[varlist]

df1$date <- as.Date(df1$date)
df1_ts <- xts(df1$Credit, df1$date)
df1_ts = as.ts(df1_ts)

# #Fitting the AR Model to the time series
# AR <- arima(df1_ts, order = c(2,0,0))
# print(AR)
# 
# ts.plot(df1_ts)
# AR_fit <- df1_ts - residuals(AR)
# points(AR_fit, type = "l", col = 2, lty = 2)

#---------
#install.packages("rucm")
modelNile <- ucm(formula = Nile~0, data = Nile, level = TRUE)
modelNile <- ucm(formula = Nile~0, data = Nile, cycle = TRUE, cycle.period = 2)

Nile1 = Nile


model1 <- ucm(formula = df1_ts~0, data = df1_ts, cycle = TRUE, cycle.period=2)

plot(df1_ts, ylab = "Flow of Nile")

lines(model1$s.level, col = "blue")

legend("topright", legend = c("Observed flow","S_level"), col = c("black","blue"), lty = 1)

plot