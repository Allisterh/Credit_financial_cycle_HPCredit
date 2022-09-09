#Merge Data
library("DataCombine")
library(dplyr)
library(vars)
#Merge Data

setwd("D:/GitHub/HPCredit/Data Collection")
sink("Granger_Correlation.txt", append=TRUE, split=TRUE) #Save output from this script


#Read raw data
df1 <- read.table("MergedData-Raw.txt", header=TRUE, sep=",")

df1 <- na.omit(df1) 


df2 <- df1 %>%
    filter(ID=="US")


# full sample
vardata0 <- df2[c("date","HHCredit_GDP_cycle_1600","HPIndex_cycle_1600")]
vardata0$date = as.Date(vardata0$date)
min(vardata0$date)

vardata0 <- df2[c("HHCredit_GDP_cycle_1600","HPIndex_cycle_1600")]
vardata0 <- ts(vardata0, start=1970, frequency=4) 

# find optimal number of lags
info.crit=VARselect(vardata0,lag.max=8,type="const")
info.crit

# estimate VAR(4)
model0=VAR(vardata0,p=5,type="const")
model1=VAR(vardata0,ic="SC",type="const")
summary(model0)
summary(model1)
# get roots (modulus)
roots(model0)
#since all values <1 -> this means the series is stationary
# model evaluation
var2.serial=serial.test(model0,lags.pt=7,type="PT.adjusted")
var2.serial
# test for serial correlation
#p values cloes to 0 is not a good news. we left something out here

# test for Granger causality
causality(model0,cause=c("HHCredit_GDP_cycle_1600"))$Granger
causality(model1,cause=c("HHCredit_GDP_cycle_1600"))$Granger


causality(model0,cause=c("HPIndex_cycle_1600"))$Granger
causality(model1,cause=c("HPIndex_cycle_1600"))$Granger

#---grangertest
vardata1 <- df[c("date","HHCredit_GDP_cycle_1600","HPIndex_cycle_1600")]
vardata1$date = as.Date(vardata1$date)
grangertest(HHCredit_GDP_cycle_1600 ~ HPIndex_cycle_1600, order = 5, data = vardata1)
grangertest(HPIndex_cycle_1600 ~ HHCredit_GDP_cycle_1600 , order = 5, data = vardata1)


#-------------lambda=400000


# full sample
vardata1 <- df2[c("date","HHCredit_GDP_cycle_400k","HPIndex_cycle_400k")]
vardata1$date = as.Date(vardata0$date)
min(vardata0$date)

vardata0 <- df2[c("HHCredit_GDP_cycle_400k","HPIndex_cycle_400k")]
vardata0 <- ts(vardata0, start=1970, frequency=4) 

# find optimal number of lags
info.crit=VARselect(vardata0,lag.max=8,type="const")
info.crit

# estimate VAR(5)
model0=VAR(vardata0,p=5,type="const")
model1=VAR(vardata0,ic="SC",type="const")
summary(model0)
summary(model1)
# get roots (modulus)
roots(model0)
#since all values <1 -> this means the series is stationary
# model evaluation
var2.serial=serial.test(model0,lags.pt=5,type="PT.adjusted")
var2.serial


var2.serial1=serial.test(model1,lags.pt=5,type="PT.adjusted")
var2.serial1
# test for serial correlation
#p values cloes to 0 is not a good news. we left something out here

# test for Granger causality
causality(model0,cause=c("HHCredit_GDP_cycle_400k"))$Granger
causality(model1,cause=c("HHCredit_GDP_cycle_400k"))$Granger


causality(model0,cause=c("HPIndex_cycle_400k"))$Granger
causality(model1,cause=c("HPIndex_cycle_400k"))$Granger


vardata1 <- df[c("date","HHCredit_GDP_cycle_400k","HPIndex_cycle_400k")]
vardata1$date = as.Date(vardata1$date)
grangertest(HHCredit_GDP_cycle_400k ~ HPIndex_cycle_400k, order = 5, data = vardata1)
grangertest(HPIndex_cycle_400k ~ HHCredit_GDP_cycle_400k , order = 5, data = vardata1)

cor(df2$HHCredit_GDP_cycle_1600,df2$HPIndex_cycle_1600)

vardata2 <- df2[c("HHCredit_GDP_cycle_400k","HPIndex_cycle_400k")]
vardata2 <- ts(vardata2, start=1970, frequency=4) 
acf(vardata2, lag.max=24, plot=FALSE)
pacf(vardata2, lag.max=24)

vardata3 <- df2[c("HHCredit_GDP_cycle_1600","HPIndex_cycle_1600")]
vardata3 <- ts(vardata3, start=1970, frequency=4) 
acf(vardata3, lag.max=24)
pacf(vardata3, lag.max=24)

dev.off()

#---Loop and correlation
library(pipeR)
library(magrittr)
df= df1

df8 <- df %>% group_by(ID) %>%
  mutate(correlation1 = cor(HHCredit_GDP_cycle_1600,HPIndex_cycle_1600)) %>%
  mutate(correlation2 = cor(HHCredit_GDP_cycle_400k,HPIndex_cycle_400k))



df9<- df8 %>% 
  group_by(borrowers_country) %>%  # Continue processing...
  #summarise(correlation1600 = mean(correlation1)) %>%
  summarise(correlation1600 = mean(correlation1),correlation400k = mean(correlation2))

print(df9, n=36)
