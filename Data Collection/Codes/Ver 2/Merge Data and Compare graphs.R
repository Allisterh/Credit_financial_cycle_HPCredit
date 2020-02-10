#Merge Data
library("DataCombine")
library(dplyr)

#Merge Data

setwd("D:/GitHub/HPCredit/Data Collection")


df2 <- read.table("HHCredit_GDP_HPfilter.txt", header=TRUE, sep=",")
df3 <- read.table("HPindex_HPfilter.txt", header=TRUE, sep=",")


df2$date=as.Date(df2$date)
df3$date=as.Date(df3$date)

#df <- merge(df1, df2, by=c("ID","date"), all=TRUE)

df <- merge(df2, df3, by=c("ID","date"), all.x=TRUE)

names(df)[4]="HHCredit"
names(df)[11]="HPIndex_trend_1600"
names(df)[12]="HPIndex_cycle_1600"
names(df)[13]="HPIndex_trend_400k"
names(df)[14]="HPIndex_cycle_400k"

write.table(df, "MergedData-Raw.txt", sep=",")


#Read raw data
df1 <- read.table("MergedData-Raw.txt", header=TRUE, sep=",")

df1 <- na.omit(df1) 



#------------------------------------------
#GRAPH 2 series Lamda = 1600

#This part of code is to shape series into one graphs
df1$date = as.Date(df1$date)

varlist1 = c("ID", "date", "borrowers_country", "HHCredit_GDP_cycle_1600", "HPIndex_GDP_cycle_1600")
df2= df1[varlist1]
df2$date = as.Date(df1$date)
names(df2)[5]="HPIndex_cycle_1600"

#head(df8)
#table(df8$variable)


#Cycle
ggplot(melt(df2, c(1,2,3)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "HPIndex and Credit",
       subtitle = "cycle - HP filter decomp | lambda=1600")
ggsave("./graphs/HPIndex_Credit_cycle_1600.pdf", width=11, height=8.5)



#GRAPH 2 series Lamda = 400000

#This part of code is to shape series into one graphs

varlist3 = c("ID", "date", "borrowers_country", "HHCredit_GDP_cycle_400k", "HPIndex_GDP_cycle_400k")
df3= df1[varlist3]
df3$date = as.Date(df3$date)
names(df3)[5]="HPIndex_cycle_400k"

#head(df8)
#table(df8$variable)


#Cycle
ggplot(melt(df3, c(1,2,3)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "HPIndex and Credit",
       subtitle = "cycle - HP filter decomp | lambda=400k")
ggsave("./graphs/HPIndex_Credit_cycle_400k.pdf", width=11, height=8.5)



