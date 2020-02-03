library("readxl")
library(reshape)
library(mFilter)
#Import HP data downloaded from Dallas Fed 

setwd("~/GitHub/HPCredit/Data Collection")
HPI = read_excel("HPI.xlsx")
PDI = read_excel("PDI.xlsx")

#View(HPI)

#Cleaning up data
HPI = HPI[-c(1),]

HPI = HPI[,-((ncol(HPI)-1):ncol(HPI))] #select last columns
#data[,(ncol(data)-2):ncol(data)]

names(HPI)[1] = "date"

HPI = as.data.frame(HPI) #melt functions need data to be as data frame
HPI = melt(HPI, id="date")

names(HPI)[3] = "HPI"

#############

View(PDI)
PDI = PDI[-c(1),]

PDI = PDI[,-((ncol(PDI)-1):ncol(PDI))] #select last columns
names(PDI)[1] = "date"

PDI = as.data.frame(PDI) #melt functions need data to be as data frame
PDI = melt(PDI, id="date")

names(PDI)[3] = "PDI"


## Merge data
HPinc <- merge(HPI, PDI, by=c("variable","date"), all.x=TRUE)

names(HPinc)[1] = "ID"

HPinc$HPinc = HPinc$HPI/HPinc$PDI

write.table(HPinc, "HPinc_full.txt", sep=",")
HPinc1=subset(HPinc, select=-c(HPI,PDI))

write.table(HPinc, "HPinc.txt", sep=",")

HPinc1 <- HPinc1%>%
  mutate(date = as.Date(as.yearqtr(date, "%Y:Q%q")))


ggplot(HPinc1, aes(date, HPinc1$HPinc, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price to Income",
       subtitle = "quarterly data")


## HP filter


HPinc1$ID = as.factor(HPinc1$ID)

df2 <- HPinc1
#df2 = as.ts(df2, c(2009, 1), end=c(2014, 12), frequency=12)

#How to make data equal length across country
library(reshape)
df3 = cast(df2, date~ID)

#clean data



#data[data[, "Var1"]>10, ]
#data[data$Var1>10, ]
#subset(data, Var1>10)
df3 = df3[df3[,"date"]>="1999-01-01",] #limit data to 1999
df3 %>% select_if(~ any(is.na(.))) %>% names()
#"CH" "CL" "CN" "ID" "IE" "IN" "LU"

#df4 = df3[,sapply(df3, function(x) !any(is.na(x)))]

#m.sapply <- function(x, ...) "attributes<-"(sapply(x, ...), attributes(x))
#m.sapply(d, function(x) x)

m.sapply <- function(x, ...) "attributes<-"(sapply(x, ...), attributes(x))
df4 = df3[,m.sapply(df3, function(x) !any(is.na(x)))]
##*remove any country with NA data


#Apply is for applying a function , sapply return a vector, lapply returns a list
write.table(df4, "HHCredit.txt", sep=",")

library(reshape)
library(mFilter)

ncol(df4) #31 countries
#remove country with NA values
#library(dplyr)
#Itun %>%
#  select_if(~ !any(is.na(.))
#Itun %>% select_if(~ any(is.na(.))) %>% names()


#str(df3)
#str(df4)

str(df4$date)
df4$date = as.character(df4$date)
df4 = as.data.frame(df4)

#Melt together
df5 = melt(df4, id.vars = "date", measure.vars = names(df4)[-1])
#df5 = melt(df3, id.vars = "date", measure.vars = names(df3)[-1])
#timeseries

#Tryout hpfilter
library(plm)

df6 <- df5 %>% group_by(variable) %>% 
  pdata.frame(., index = c("variable","date")) %>% 
  mutate(HHCredit_GDP_gap = value-mFilter::hpfilter(value, type = "lambda", freq = 400000)$trend)

names(df6)[2] = "ID"
names(df6)[4] = "HP_inc_gap"

write.table
write.table(df6[,c(1,2,4)], "HP_Inc_gap.txt", sep=',' )

df6$date = as.Date(df6$date)

#Graphing
ggplot(df6, aes(date, HP_inc_gap, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Price",
       subtitle = "as percentage of disposible income - HP decomp cycle | lambda=400000")

#Change name to ID
df7 = df6
df7 = as.data.frame(df7)
names(df7)[2]="country"
df7$country0 = as.character(df7$country)


from=c("Australia", "Belgium", "Canada" , "Switzerland",  "Germany", "Denmark"  ,   "Spain"    ,   "Finland"   , "France" ,     "UK"       , "Ireland"  ,   "Italy"   ,    "Japan"     ,  "S. Korea" ,   "Luxembourg" , "Netherlands", "Norway"    ,  "New Zealand" ,"Sweden"     , "US"     , "S. Africa" ,  "Croatia"    , "Israel")
to=c("AU", "BE", "CA" , "CH",  "DE", "DK"  ,   "ES"    ,   "FI"   , "FR" ,     "GB"       , "IE"  ,   "IT"   ,    "JP"     ,  "KR" ,   "LU" , "NL", "NO"    ,  "NZ" ,"SE"     , "US"     , "ZA" ,  "HR"    , "IR")

map = setNames(to,from)
df7$country0[] = map[df7$country0]
View(df7)

names(df7)[ncol(df7)]="ID"

rownames(df7) <- c()

write.table(df7[c("date","ID","HP_inc_gap")], "HP_inc_gap.txt", sep = ',',row.names = FALSE)
#Work on later
#Disposible Income
#https://stats.oecd.org/viewhtml.aspx?datasetcode=NAAG&lang=en#
#https://stats.oecd.org/restsdmx/sdmx.ashx/GetData/NAAG/AUS+AUT+BEL+CAN+CHL+CZE+DNK+EST+FIN+FRA+DEU+GRC+HUN+ISL+IRL+ISR+ITA+JPN+KOR+LVA+LTU+LUX+MEX+NLD+NZL+NOR+POL+PRT+SVK+SVN+ESP+SWE+CHE+TUR+GBR+USA+EMU+OTO+NMEC+BRA+CHN+COL+CRI+IND+IDN+RUS+ZAF.GDPCPC+B6GS14_S15HCPC+B7GS14_S15HCPC/all?startTime=2000&endTime=2018
