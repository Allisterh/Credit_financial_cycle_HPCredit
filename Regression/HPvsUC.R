 rm(list=ls())

#Merge Data
library("DataCombine")
library(dplyr)
library(reshape2)
library(ggplot2)
library(fredr)
 library(extrafont)
#Merge Data

#Version selections#####
ver='VAR_2_crosscycle'
country = 'GB'

working_dir=sprintf("D:/GitHub/HPCredit/Regression/%s/R",ver)
setwd(working_dir) 

#automate file name  paste(x, ".mean", sep="")

#Credit and HPI merge

#Read raw data
HP_filepath = sprintf("../../../Data Collection/1.Latest/HPindex_HPfilter_%s.txt",country)
df1 <- read.table(HP_filepath, header=TRUE, sep=",")

# df1 <- na.omit(df1[-c(2)]) Remove country column

Credit_filepath = sprintf("../../../Data Collection/1.Latest/Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy

df <- merge(df1, df2, by=c("ID","date"))
df <-subset(df, date>as.Date("1988-12-30"))
df <-subset(df, date<as.Date("2020-01-01"))


## Recession

add_rec_shade<-function(st_date,ed_date,shade_color="darkgray")
{
  library(fredr)
  library(ecm)
  library(ggplot2)
  fredr_set_key("809225ae418303bc5fd2979e182aa537")
  
  st_date<-as.Date("1988-12-30")
  ed_date<-as.Date("2020-01-01")
  
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

df3 <- read.table(sprintf("../Output/OutputData/uc_yc_%s.txt",country), header=FALSE, sep=",")
df = cbind(df,df3)
df$date = as.Date(df$date) 

#head(df8)
#table(df8$variable)

#Cycles var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "HPIndex_HPcycle", "V2", "V4")

#Credit Cycle var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "V2")

df6 = df[varlist2]
names(df6)[4]="UC_Credit_Cycle"
names(df6)[3]="HP_Credit_Cycle"

p1<-ggplot(melt(df6, c(1,2)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Credit cycle: %s",country)  )
ggsave( sprintf("../Output/graphs/Credit_cycle_%s.pdf",country) , width=8, height=5)

#HP Cycle
varlist3 = c("ID", "date", "HPIndex_HPcycle", "V4")
df7 = df[varlist3]
names(df7)[4]="UC_HPI_Cycle"
names(df7)[3]="HP_HPI_Cycle"


p2<- ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(x = NULL, y = NULL,
       title = sprintf("Housing Price cycle: %s",country)  )
ggsave( sprintf("../Output/graphs/HP_cycle_%s.pdf",country) , width=8, height=5)

#Trends 
#Credit Trends

varlist4 = c("ID", "date", "Credit_HPtrend", "Credit_log", "V1")
df7 = df[varlist4]
names(df7)[3]="HP_Credit_Trend"
names(df7)[4]="Series"
names(df7)[5]="UC_Credit_Trend"
df7[c(3:5)] = exp(df7[c(3:5)]/100)

p3<-ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Credit Trend: %s , as percentage of GDP",country)  )
ggsave( sprintf("../Output/graphs/Credit_Trend_%s.pdf",country) , width=8, height=5)


#Housing Price Index Trends

varlist5 = c("ID", "date", "HPIndex_HPtrend", "HPIndex_log", "V3")
df7 = df[varlist5]
names(df7)[3]="HP_Credit_Trend"
names(df7)[4]="Series"
names(df7)[5]="UC_Credit_Trend"
df7[c(3:5)] = exp(df7[c(3:5)]/100)

p4<-ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  add_rec_shade(min(fred_data$date),max(fred_data$date))+
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme(legend.position = "bottom") +
  theme_light() +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(x = NULL, y = NULL,
    title = sprintf("Housing Price Index Trend: %s , Index 2010=100",country)  )
ggsave( sprintf("../Output/graphs/HP_Trend_%s.pdf",country) , width=8, height=5)

library(patchwork)
(p1|p3)/(p2|p4)
ggsave( sprintf("../Output/graphs/HP_Credit_4graphs_%s.pdf",country) , width=8, height=5)
