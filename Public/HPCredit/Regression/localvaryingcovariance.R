rm(list=ls())

#Import data
country = 'US'

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

VAR2x1_filepath = sprintf("Bayesian_UC_VAR2_drift_Crosscycle1lag/OutputData/filter_uc_%s.csv",country)

df <- read.table(VAR2x1_filepath, header=FALSE, sep=",")
cor(df[,7],df[,8])


country = 'UK'

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

VAR2x1_filepath = sprintf("Bayesian_UC_VAR2_drift_Crosscycle1lag/OutputData/filter_uc_%s.csv",country)

df <- read.table(VAR2x1_filepath, header=FALSE, sep=",")
cor(df[,7],df[,8])