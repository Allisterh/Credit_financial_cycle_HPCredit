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
ver='VAR_2'
country = 'US'

working_dir=sprintf("D:/GitHub/HPCredit/Regression/%s/R",ver)
setwd(working_dir) 

#automate file name  paste(x, ".mean", sep="")

#Credit and HPI merge

#Read raw data
HP_filepath = sprintf("../../../Data Collection/1.Latest/Paper2/MergedData_Matlab_%s.txt",country)
df1 <- read.table(HP_filepath, header=TRUE, sep=",")

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


df4 <- read.table(sprintf("../../AR_2/Output/OutputData/uc_yc_credit_%s.txt",country), header=FALSE, sep=",")
names(df4)[2] = "UC_Univariate_Credit"
df5 <- read.table(sprintf("../../AR_2/Output/OutputData/uc_yc_hpi_%s.txt",country), header=FALSE, sep=",")
names(df5)[2] = "UC_Univariate_HPI"

df = cbind(df,df4[2],df5[2])

#head(df8)
#table(df8$variable)

#Cycles var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "HPIndex_HPcycle", "V2", "V4")

#Credit Cycle var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "V2","UC_Univariate_Credit")

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
ggsave( sprintf("../../AR_2/Output/graphs/Credit_cycle_%s.pdf",country) , width=8, height=5)

p1


#HPI Cycle

varlist3 = c("ID", "date", "HPIndex_HPcycle", "V4", "UC_Univariate_HPI")
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
ggsave( sprintf("../../AR_2/Output/graphs/HP_cycle_%s.pdf",country) , width=8, height=5)

p2
library(patchwork)
(p1)/(p2)
ggsave( sprintf("../../AR_2/Output/graphs/HP_Credit_2graphs_%s.pdf",country) , width=8, height=8)
