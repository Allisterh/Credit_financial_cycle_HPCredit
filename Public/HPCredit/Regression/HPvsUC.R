rm(list=ls())

#Merge Data
# library("DataCombine")
library(dplyr)
library(reshape2)
library(ggplot2)
library(fredr)
#library(extrafont)
#Merge Data

#Version selections#####
ver='Bayesian_UC_VAR2_drift_Crosscycle2lags'
country = 'US'

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
getwd()

#automate file name  paste(x, ".mean", sep="")

#Credit and HPI merge

#Read raw data
# HP_filepath = sprintf("../../../Data Collection/1.Latest/Paper1/creditHPI_Matlab_US_%s.csv",country)
# df1 <- read.table(HP_filepath, header=TRUE, sep=",")

# # df1 <- na.omit(df1[-c(2)]) Remove country column

# Credit_filepath = sprintf("../../../Data Collection/1.Latest/Paper1/Credit_HPfilter_%s.csv",country)
# df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
# df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy
filepath = sprintf("../Data Collection/1.Latest/Paper1/creditHPI_Matlab_%s.csv",country)
df <- read.table(filepath, header=FALSE, sep=",")
df <- df[-c(1:2),]
df <- df[,-c(3:4)]
df[,3] <- exp(df[,3]/100)
df[,5] <- exp(df[,5]/100)

st_date<-as.Date("1991-04-01")
ed_date<-as.Date("2021-07-01")
V <- seq(as.Date(st_date), as.Date(ed_date), by="quarters")
df$date <- V

# df <- merge(df1, df2, by=c("ID","date"))
# df <-subset(df, date>=as.Date("1991-04-01"))
# df <-subset(df, date<=as.Date("2021-07-01"))


## Recession

add_rec_shade<-function(st_date,ed_date,shade_color="darkgray")
{
  library(fredr)
  library(ecm)
  library(ggplot2)
  fredr_set_key("809225ae418303bc5fd2979e182aa537")
  
  st_date<-as.Date("1988-12-30")
  ed_date<-as.Date("2021-07-01")
  
  ## US: USRECD
  ## UK: GBRRECD
  recession<-fredr(series_id = "USRECD",observation_start = as.Date(st_date),observation_end = as.Date(ed_date))
  
  recession$diff<-recession$value-lagpad(recession$value,k=1)
  recession<-recession[!is.na(recession$diff),]
  recession.start<-recession[recession$diff==1,]$date
  recession.end<-recession[recession$diff==(-1),]$date
  
  if(length(recession.start)>length(recession.end))
  {recession.end<-c(recession.end,Sys.Date())}
  if(length(recession.end)>length(recession.start))
  {recession.start<-c(min(recession$date),recession.start)}
  
  recs<-as.data.frame(cbind(recession.start,recession.end))
  recs$recession.start<-as.Date(as.numeric(recs$recession.start),origin=as.Date("1970-01-01"))
  recs$recession.end<-as.Date(recs$recession.end,origin=as.Date("1970-01-01"))
  if(nrow(recs)>0)
  {
    rec_shade<-geom_rect(data=recs, inherit.aes=F, 
                         aes(xmin=recession.start, xmax=recession.end, ymin=-Inf, ymax=+Inf), 
                         fill=shade_color, alpha=0.5)
    return(rec_shade)
  }
}
#------------------------------------------
#GRAPH 2 series Lamda = 1600

#This part of code is to shape series into one graphs
# df1$date = as.Date(df1$date)
# 
# varlist1 = c("ID", "date", "borrowers_country", "HHCredit_GDP_cycle_1600", "HPIndex_GDP_cycle_1600", "HHCredit_GDP_trend_1600", "HPIndex_GDP_trend_1600")
# df2= df1[varlist1]
# df2$date = as.Date(df2$date)
# names(df2)[5]="HPIndex_cycle_1600"
# names(df2)[7]="HPIndex_trend_1600"
# 
# df3 <- df2 %>%
#   filter(ID=="US")
# df3 <-na.omit(df3)

#Merge data with UC model data

filepath = sprintf("%s/OutputData/",ver)
filepath2 = sprintf("filter_uc_%s.csv",country)
filepath <- paste(filepath,filepath2, sep = "")

df3 <- read.table(filepath, header=FALSE, sep=",")
df3 <- df3[,c(1:2,4:5,7:8)]
df3[,1] <- exp(df3[,1]/100)
df3[,3] <- exp(df3[,3]/100)

df = cbind(df,df3)
names(df) <- c("Credit series", "HPI series", "Credit_Trend_HP", "Credit_Cycle_HP",
   "HPI_Trend_HP", "HPI_Cycle_HP", "date", "Credit_Trend_UC", "Credit_Cycle_UC",
   "HPI_Trend_UC", "HPI_Cycle_UC", "Credit local growth rate", "HPI local growth rate")
#df$date = as.Date(df$date) 

#head(df8)
#table(df8$variable)

#Cycles var name list
varlist = c("date", "Credit_Cycle_HP", "HPI_Cycle_HP", "Credit_Cycle_UC", "HPI_Cycle_UC")

#Credit Cycle var name list
varlist = c("date", "Credit_Cycle_HP", "Credit_Cycle_UC")

df6 = df[varlist]

p1<-ggplot(melt(df6, c(1)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Credit cycle: %s",country))


#HP Cycle
varlist = c("date", "HPI_Cycle_HP", "HPI_Cycle_UC")
df7 = df[varlist]

p2<- ggplot(melt(df7, c(1)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(x = NULL, y = NULL,
       title = sprintf("Housing Price cycle: %s",country)  )
#ggsave( sprintf("../Output/graphs/HPI_Cycle_%s.pdf",country) , width=8, height=5)

#Trends 
#Credit Trends

varlist = c("date", "Credit_Trend_HP", "Credit_Trend_UC", "Credit series")

df7 = df[varlist]

p3<-ggplot(melt(df7, c(1)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Credit Trend: %s , as percentage of GDP",country))
#ggsave( sprintf("../Output/graphs/Credit_Trend_%s.pdf",country) , width=8, height=5)


#Housing Price Index Trends
varlist = c("date", "HPI_Trend_HP", "HPI_Trend_UC", "HPI series")
df7 = df[varlist]

p4<-ggplot(melt(df7, c(1)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme(legend.position = "bottom") +
  theme_light() +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(x = NULL, y = NULL,
    title = sprintf("Housing Price Index Trend: %s , Index 2010=100",country))

filepath = sprintf("%s/OutputData/graphs/",ver)
filepath2 = sprintf("HP_Credit_4graphs_%s.pdf",country)
filepath <- paste(filepath,filepath2, sep = "")
library(patchwork)
(p1|p3)/(p2|p4) 
ggsave( filepath , width=8, height=5)
