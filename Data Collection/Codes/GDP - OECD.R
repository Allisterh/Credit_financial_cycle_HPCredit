library(OECD)
library(IMFData)
library(imfr)
#library(fredr)
#library(plyr)
library(dplyr)
library(quantmod)
library(ggplot2)

setwd("D:/OD/OneDrive/Credit/Laptop/HPCredit/Data Collection")

#Set dates
startdate='1969-01-01'
enddate='2019-12-01'

#Guide here https://github.com/expersso/OECD
#Search OECD
get_datasets <- OECD::get_datasets

#GDP
dataset_list <- get_datasets()
x = search_dataset("Quarterly national account", data = dataset_list)
#View(x)

dstruc <- get_data_structure("QNA")
str(dstruc, max.level = 1)
dstruc$VAR_DESC
#View(dstruc$MEASURE)

#GDP (VOBARSA: National currency, volume estimates, OECD reference)
df1 <- get_dataset("QNA",
                   filter = "AUT+AUS+BEL+BRA+CAN+CHE+CHL+CHN+COL+CRI+CZE+DEU+DNK+ESP+EST+FIN+FRA+GBR+GRC+HUN+IND+IDN+IRL+ISL+ISR+ITA+JPN+KOR+LTU+LUX+LVA+MEX+NLD+NOR+NZL+POL+PRT+ROU+RUS+SAU+SRB+SVK+SVN+SWE+TUR+USA+ZAF.B1_GE.HCPCARSA",
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

df <- df %>%
  mutate(GDP = log(GDP))

from=c("AUT","AUS","BEL","BRA","CAN","CHE","CHL","CHN","COL","CRI","CZE","DEU","DNK","ESP","EST","FIN","FRA","GBR","GRC","HUN","IND","IDN","IRL","ISL","ISR","ITA","JPN","KOR","LTU","LUX","LVA","MEX","NLD","NOR","NZL","POL","PRT","ROU","RUS","SAU","SRB","SVK","SVN","SWE","TUR","USA","ZAF")
to=c("AT","AU","BE","BR","CA","CH","CL","CN","CO","CR","CZ","DE","DK","ES","EE","FI","FR","GB","GR","HU","ID","IN","IE","IS","IL","IT","JP","KR","LT","LU","LV","MX","NL","NO","NZ","PL","PT","RO","RU","SA","RS","SK","SI","SE","TR","US","ZA")
map = setNames(to,from)
df$ID[] = map[df$ID]

str(from)
str(to)

#\r\n to replace new line with comma

write.table(df, "/Data/GDP.txt", sep=",")


##Graphing
ggplot(df, aes(date, GDP, color = ID)) +
  geom_hline(yintercept = 8, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "log GDP percapita",
       subtitle = "quarterly data")