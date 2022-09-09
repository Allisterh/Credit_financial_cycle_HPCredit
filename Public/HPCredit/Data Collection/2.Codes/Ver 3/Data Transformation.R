rm(list=ls())

#Merge Data
#library("DataCombine")
#library(quantmod)
#library(stats)
#library("tidyr")
library(dplyr)
library(reshape2)


##1. Merge Data
country = 'US'

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
setwd("../../1.Latest/Paper1")

Credit_filepath = sprintf("Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy
 

HP_filepath = sprintf("HPindex_HPfilter_%s.txt",country)
df1 <- read.table(HP_filepath, header=TRUE, sep=",")


df <- merge(df2, df1, by=c("ID","date"))
df <- subset(df, date>as.Date("1985-06-30"))
df <- subset(df, date<as.Date("2020-01-01"))

filepath = sprintf("MergedData_%s.txt",country)
write.table(df, filepath, sep=",")

varlist = c("Credit", "HPIndex", "Credit_HPcycle","HPIndex_HPcycle","Credit_HPtrend", "HPIndex_HPtrend")
df_matlab <-df[varlist]
filepath = sprintf("MergedData_Matlab_%s.txt",country)
write.table(df_matlab, filepath, sep=",")


# mean(df$HPIndex_log)
# max(df$HPIndex_log)
# min(df$HPIndex_log)
# 
# mean(df$Credit_log)
# max(df$Credit_log)
# min(df$Credit_log)
# 
# cor(df$HPIndex_log,df$Credit_log)

## Pseudo-code
##Transform into diff
df$date <- as.Date(df$date)
df_xts <- as.xts(df$Credit, df$date)
names(df_xts)<-"Credit"
df_xts$Credit_growthrate <- (df_xts$Credit / stats::lag(df_xts$Credit, 1) - 1)*100
##Summarize variable.names
df_xts <- na.omit(df_xts)
mean(df_xts$Credit_growthrate)
max(df_xts$Credit_growthrate)
min(df_xts$Credit_growthrate)


## HPIndex
df_xts <- xts(df$HPIndex, df$date)
names(df_xts)<-"HPIndex"
df_xts$HPIndex_growthrate <- (df_xts$HPIndex / stats::lag(df_xts$HPIndex, 1) - 1)*100
##Summarize variable.names
df_xts <- na.omit(df_xts)
mean(df_xts$HPIndex_growthrate)
max(df_xts$HPIndex_growthrate)
min(df_xts$HPIndex_growthrate)

##Correlation
## create 4 lag series
df_xts <- as.xts(df, order.by=df$date)
df_xts$Credit_log_1 <- stats::lag(df_xts$Credit_log, 1)
df_xts$Credit_log_2 <- stats::lag(df_xts$Credit_log, 2)
df_xts$HPIndex_log_1 <- stats::lag(df_xts$HPIndex_log, 1)
df_xts$HPIndex_log_2 <- stats::lag(df_xts$HPIndex_log, 2)

## 9 corr statistics
df_xts <- as.data.frame(df_xts)
df_xts <- na.omit(df_xts)
df_xts <- df_xts[c("Credit_log", "Credit_log_1", "Credit_log_2", "HPIndex_log", "HPIndex_log_1", "HPIndex_log_2")]

df_xts$Credit_log <- as.numeric(df_xts$Credit_log)
df_xts$Credit_log_1 <- as.numeric(df_xts$Credit_log_1)
df_xts$Credit_log_2 <- as.numeric(df_xts$Credit_log_2)
df_xts$HPIndex_log <- as.numeric(df_xts$HPIndex_log)
df_xts$HPIndex_log_1 <- as.numeric(df_xts$HPIndex_log_1)
df_xts$HPIndex_log_2 <- as.numeric(df_xts$HPIndex_log_2)

cor(df_xts$Credit_log,df_xts$Credit_log_1)
cor(df_xts$Credit_log,df_xts$Credit_log_2)
cor(df_xts$Credit_log,df_xts$HPIndex_log_1)
cor(df_xts$Credit_log,df_xts$HPIndex_log_2)


cor(df_xts$HPIndex_log,df_xts$Credit_log_1)
cor(df_xts$HPIndex_log,df_xts$Credit_log_2)
cor(df_xts$HPIndex_log,df_xts$HPIndex_log_1)
cor(df_xts$HPIndex_log,df_xts$HPIndex_log_2)

cor(df_xts$HPIndex_log,df_xts$Credit_log)