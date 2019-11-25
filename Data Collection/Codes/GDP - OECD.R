library(OECD)
library(IMFData)
library(imfr)
#library(fredr)
#library(plyr)
library(dplyr)
library(quantmod)
library(ggplot2)

setwd("D:/OD/OneDrive/Project/HPCredit/Data Collection")

#Set dates
startdate='1969-01-01'
enddate='2019-12-01'

#Guide here https://github.com/expersso/OECD
#Search OECD
get_datasets <- OECD::get_datasets

#GDP
#dataset_list <- get_datasets()
#x = search_dataset("Quarterly national account", data = dataset_list)
#View(x)

dstruc <- get_data_structure("QNA")
str(dstruc, max.level = 1)
dstruc$VAR_DESC
dstruc$LOCATION
#View(dstruc$MEASURE)
#Export table of QNA measures
write.table(dstruc$MEASURE, "QNAMeasures.txt", sep=",")
write.table(dstruc$LOCATION, "QNALocations.txt", sep=",")

#List of available call codes for GDP: in million unit

#"5","CQR","National currency, current prices, quarterly levels"
#"6","CQRSA","National currency, current prices, quarterly levels, seasonally adjusted"
#"8","IND","Volume and price indices"
#"10","GPSA","Growth rate compared to previous quarter, seasonally adjusted"
#"11","GRW","Growth rates"
#"12","GYSA","Growth rate compared to the same quarter of previous year, seasonally adjusted"
#"18","LNBARSA","National currency, chained volume estimates, national reference year, annual levels, seasonally adjusted"
#"19","VOL","Volumes"
#"20","LNBQR","National currency, chained volume estimates, national reference year, quarterly levels"
#"21","LNBQRSA","National currency, chained volume estimates, national reference year, quarterly levels, seasonally adjusted"
#"22","VIXNB","Volume index, national base/reference year"
#"23","VIXNBSA","Volume index, national base/reference year, seasonally adjusted"
#"24","VIXOBSA","Volume index, OECD reference year, seasonally adjusted"
#"25","VNBAR","National currency, constant prices, national base year, annual levels"
#"26","VNBARSA","National currency, constant prices, national base year, annual levels, seasonally adjusted"
#"27","VNBQR","National currency, constant prices, national base year, quarterly levels"
#"28","VNBQRSA","National currency, constant prices, national base year, quarterly levels, seasonally adjusted"
#"29","VOBARSA","National currency, volume estimates, OECD reference year, annual levels, seasonally adjusted"
#"30","VPVOBARSA","US dollars, volume estimates, fixed PPPs, OECD reference year, annual levels, seasonally adjusted"
#"31","PERSA","Persons, seasonally adjusted"
#"32","PER","Persons"
#"34","HCPCARSA","Per Head, US $, current prices, current PPPs, seasonally adjusted"
#"35","HVPVOBARSA","Per Head, US $, constant prices, fixed PPPs, OECD reference year, seasonally adjusted"

#Potential pick for this exercise: Nominal GDP
#"5","CQR","National currency, current prices, quarterly levels"



#Query is too long, will need to shorten it
#GDP (VOBARSA: National currency, volume estimates, OECD reference)
df1 <- get_dataset("QNA",
                   filter = "AUT+AUS+BEL+BRA+CAN+CHE+CHL+CHN+COL+CRI+CZE+DEU+DNK+ESP+EST+FIN+FRA+GBR+GRC+HUN+IND+IDN+IRL+ISL+ISR+ITA+JPN+KOR+LTU+LUX+LVA+MEX+NLD+NOR+NZL+POL+PRT+ROU+RUS+SAU+SVK+SVN+SWE+TUR+USA+ZAF.B1_GE.CQR",
                   start_time = startdate, end_time = enddate,
                   pre_formatted = TRUE)

df1 <- df1%>%
  filter(grepl("^(Q)", FREQUENCY))%>%
  mutate(date = as.Date(as.yearqtr(obsTime, "%Y-Q%q")))
  
#Saving data

myvars <- c("LOCATION", "date", "obsValue")
df <- df1[myvars]


names(df)[1]<-"ID"
names(df)[3]<-"GDP"

#df <- df %>%
#  mutate(GDP = log(GDP))

from=c("AUT","AUS","BEL","BRA","CAN","CHE","CHL","CHN","COL","CRI","CZE","DEU","DNK","ESP","EST","FIN","FRA","GBR","GRC","HUN","IND","IDN","IRL","ISL","ISR","ITA","JPN","KOR","LTU","LUX","LVA","MEX","NLD","NOR","NZL","POL","PRT","ROU","RUS","SAU","SRB","SVK","SVN","SWE","TUR","USA","ZAF")
to=c("AT","AU","BE","BR","CA","CH","CL","CN","CO","CR","CZ","DE","DK","ES","EE","FI","FR","GB","GR","HU","ID","IN","IE","IS","IL","IT","JP","KR","LT","LU","LV","MX","NL","NO","NZ","PL","PT","RO","RU","SA","RS","SK","SI","SE","TR","US","ZA")
map = setNames(to,from)
df$ID[] = map[df$ID]

str(from)
str(to)

#\r\n to replace new line with comma

write.table(df, "GDP.txt", sep=",")
table(df$ID)
length(table(df$ID))

##Graphing
ggplot(df, aes(date, GDP, color = ID)) +
  geom_hline(yintercept = 8, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Nominal GDP",
       subtitle = "quarterly data")


#Log graphing
df2 <- df %>%
  mutate(GDP = log(GDP))

ggplot(df2, aes(date, GDP, color = ID)) +
  geom_hline(yintercept = 8, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Log Nominal GDP",
       subtitle = "quarterly data - OECD")
