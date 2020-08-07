#Merge Data
library("DataCombine")
library(dplyr)
library(reshape2)
library(ggplot2)
#Merge Data

#Version selections
ver='VAR_2'
country = 'US'

working_dir=sprintf("D:/GitHub/HPCredit/Regression/%s/R",ver)
setwd(working_dir) 

#automate file name  paste(x, ".mean", sep="")

#Credit and HPI merge

#Read raw data
HP_filepath = sprintf("../../../Data/Input/HPindex_HPfilter_%s.txt",country)
df1 <- read.table(HP_filepath, header=TRUE, sep=",")

# df1 <- na.omit(df1[-c(2)]) Remove country column

Credit_filepath = sprintf("../../../Data/Input/Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy

df <- merge(df1, df2, by=c("ID","date"))

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

df = df[-nrow(df),]
df = df[-1,]
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

ggplot(melt(df6, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
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


ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Housing Price cycle: %s",country)  )
ggsave( sprintf("../Output/graphs/HP_cycle_%s.pdf",country) , width=8, height=5)

#Trends 
#Credit Trends

varlist4 = c("ID", "date", "Credit_HPtrend", "V1", "Credit_log")
df7 = df[varlist4]
names(df7)[4]="UC_Credit_Trend"
names(df7)[3]="HP_Credit_Trend"
names(df7)[5]="Series"
df7[c(3:5)] = exp(df7[c(3:5)]/100)

ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Credit Trend: %s , as percentage of GDP",country)  )
ggsave( sprintf("../Output/graphs/Credit_Trend_%s.pdf",country) , width=8, height=5)


#Housing Price Index Trends

varlist5 = c("ID", "date", "HPIndex_HPtrend", "V3", "HPIndex_log")
df7 = df[varlist5]
names(df7)[4]="UC_Credit_Trend"
names(df7)[3]="HP_Credit_Trend"
names(df7)[5]="Series"
df7[c(3:5)] = exp(df7[c(3:5)]/100)

ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
    title = sprintf("Housing Price Index Trend: %s , Index 2010=100",country)  )
ggsave( sprintf("../Output/graphs/HP_Trend_%s.pdf",country) , width=8, height=5)


