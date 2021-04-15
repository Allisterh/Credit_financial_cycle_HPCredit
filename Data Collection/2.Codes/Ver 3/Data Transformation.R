rm(list=ls())

#Merge Data
library("DataCombine")
library(dplyr)
library(reshape2)

##1. Merge Data
country = 'US'
setwd("D:/GitHub/HPCredit/Data Collection/1.Latest")

Credit_filepath = sprintf("Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy


HP_filepath = sprintf("HPindex_HPfilter_%s.txt",country)
df1 <- read.table(HP_filepath, header=TRUE, sep=",")


df <- merge(df2, df1, by=c("ID","date"))
df <-subset(df, date>as.Date("1988-06-30"))
df <-subset(df, date<as.Date("2020-01-01"))

filepath = sprintf("MergedData_%s.txt",country)
write.table(df, filepath, sep=",")

varlist = c("Credit", "HPIndex", "Credit_HPcycle","HPIndex_HPcycle","Credit_HPtrend", "HPIndex_HPtrend")
df_matlab <-df[varlist]
filepath = sprintf("MergedData_Matlab_%s.txt",country)
write.table(df_matlab, filepath, sep=",")


mean(df$HPIndex_log)
max(df$HPIndex_log)
min(df$HPIndex_log)

mean(df$Credit_log)
max(df$Credit_log)
min(df$Credit_log)

cor(df$HPIndex_log,df$Credit_log)
