library(dplyr)
library(ggplot2)
library(zoo)
library(BIS)

setwd("~/GitHub/HPCredit/Data Collection")

datasets <- BIS::get_datasets()
#head(datasets, 20)

#country list (with HP available) 
clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|SA|RS|SK|SI|SE|TR|US|ZA
)"

rates <- get_bis(datasets$url[datasets$name == "Credit to the non-financial sector"], quiet = TRUE)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl(clist, borrowers_cty))%>%
  filter(grepl("^(H)", tc_borrowers))%>%
  filter(grepl("^(All sectors)", lending_sector))%>%
  filter(grepl("^(770)", unit_type))

#%>%
#filter(grepl("^(U)", tc_adjust)) 

#%>%
#group_by(borrowers_cty) %>%
#mutate(growth = c(NA,diff(obs_value, lag = 4))*100/obs_value)

table(rates_plot$ref_area)
table(rates_plot$date)
table(rates_plot$unit_type)
table(rates_plot$tc_adjust)

ggplot(rates_plot, aes(date, obs_value, color = borrowers_country)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Credit to household",
       subtitle = "as percentage of GDP")


#Saving Data

myvars <- c("borrowers_cty", "date", "obs_value")
df <- rates_plot[myvars]

names(df)[1]<-"ID"
names(df)[3]<-"HHCredit"

write.table(df, "HHCredit.txt", sep=",")

df$ID = as.factor(df$ID)

df2 <- df
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
  mutate(HHCredit_GDP_gap = mFilter::hpfilter(value, type = "lambda", freq = 3000)$cycle)

###############Non financial credit

rates <- get_bis(datasets$url[datasets$name == "Credit to the non-financial sector"], quiet = TRUE)


rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl(clist, borrowers_cty))%>%
  filter(grepl("^(P)", tc_borrowers))%>%
  filter(grepl("^(All sectors)", lending_sector))%>%
  filter(grepl("^(770)", unit_type))

#%>%
#filter(grepl("^(U)", tc_adjust)) 

#%>%
#group_by(borrowers_cty) %>%
#mutate(growth = c(NA,diff(obs_value, lag = 4))*100/obs_value)

table(rates_plot$ref_area)
table(rates_plot$date)
table(rates_plot$unit_type)
table(rates_plot$tc_adjust)

ggplot(rates_plot, aes(date, obs_value, color = borrowers_country)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Credit to nonfinancial sector",
       subtitle = "as percentage of GDP")


#Saving Data

myvars <- c("borrowers_cty", "date", "obs_value")
df <- rates_plot[myvars]

names(df)[1]<-"ID"
names(df)[3]<-"NFECredit"

write.table(df, "NFECredit.txt", sep=",")


df$ID = as.factor(df$ID)

df2 <- df
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
  mutate(HHCredit_GDP_gap = mFilter::hpfilter(value, type = "lambda", freq = 3000)$cycle)


names(df6)[2] = "ID"

write.table
write.table(df6[,c(1,2,4)], "NFECredit_GDP_gap.txt", sep=',' )