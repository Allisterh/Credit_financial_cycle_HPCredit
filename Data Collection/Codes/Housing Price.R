library(dplyr)
library(ggplot2)
library(zoo)
library(BIS)

setwd("D:/GitHub/HPCredit/Data Collection")

#Preparing Dataset
detach("package:OECD", unload = TRUE)
get_datasets = BIS::get_datasets
datasets <- get_datasets()
head(datasets, 20)

#country list
#clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|SA|RS|SK|SI|SE|TR|US|ZA)"

#modified clist to match HP from BIS
clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|FI|FR|GB|GR|HU|IE|ID|IN|IL|IT|JP|KR|LU|MX|NL|NO|NZ|PL|PT|RU|SA|SE|TR|US)"


#-----------------------
#House Price
#-----------------------

rates <- get_bis(datasets$url[datasets$name == "Property prices: selected series"], quiet = TRUE)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl(clist, ref_area))%>%
  filter(grepl("^(R)", value))%>%
  filter(grepl("^(Index, 2010 = 100)", unit_of_measure)) 

#771 vs 628

#table(rates_plot$ref_area)
#table(rates_plot$date)
#table(rates_plot$freq)
#table(rates$value)
#table(rates_plot$unit_of_measure)

table(rates$unit_of_measure)

ggplot(rates_plot, aes(date, obs_value, color = reference_area)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Real House Price Growth Rate",
       subtitle = "quarterly data")

write.table(rates_plot, "HP.txt", sep=",")
table(df$ID)
length(table(df$ID))

#as.yearqtr("2001 Q2")
#as.yearqtr("2001-q2","%Y-q%q") # same


#-----------------------
#Credit to household - Credit to the non-financial sector
#-----------------------

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

ggplot(rates_plot, aes(date, growth, color = borrowers_country)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Credit to household",
       subtitle = "Year to Year % change")

# Test for US Data
rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl("^(US)", borrowers_cty))%>%
  filter(grepl("^(H)", tc_borrowers))%>%
  filter(grepl("^(All sectors)", lending_sector))%>%
  filter(grepl("^(XDC)", unit_type))%>%
  filter(grepl("^(U)", tc_adjust))%>%
  group_by(borrowers_cty) %>%
  mutate(growth = c(NA,diff(obs_value, lag = 4))*100/obs_value)

table(rates_plot$ref_area)
table(rates_plot$date)
table(rates_plot$unit_type)
table(rates_plot$tc_adjust)

ggplot(rates_plot, aes(date, growth, color = borrowers_country)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Credit to household",
       subtitle = "Year to Year % change")

#test completed


library(tidyverse)
library(AER)
data("Grunfeld",package = "AER")
Grunfeld1 <- Grunfeld %>% 
  select(firm,year,invest) %>%
  group_by(firm) %>%
  mutate(growth = c(NA,diff(invest, lag = 4))/invest)

Grunfeld2 <- Grunfeld %>% 
  select(firm,year,invest) %>%
  group_by(firm) %>%
  mutate(growth = c(NA,diff(invest))/invest)


View(Grunfeld1)
View(Grunfeld2)


#Full sample HP filter

myvars <- c("ref_area", "date", "obs_value")
df <- rates_plot[myvars]

names(df)[1]<-"ID"
names(df)[3]<-"HHCredit"

names(df)[1] = "variable"
names(df)[3] = "value"

df7 <- df %>% group_by(variable) %>% 
  pdata.frame(., index = c("variable","date")) %>%
  mutate(HHCredit_GDP_trend = mFilter::hpfilter(value, type = "lambda", freq = 400000)$trend)%>%
  mutate(HHCredit_GDP_gap = value-HHCredit_GDP_trend)

names(df7)[1] = "ID"
names(df7)[2] = "date"
names(df7)[4] = "HP_trend"
names(df7)[5] = "HP_gap"

df7$date = as.Date(df7$date)

df7$date = as.Date(df7$date)
ggplot(df7, aes(date, HP_gap, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price index 2010 = 100",
       subtitle = "Detrended Gap - HPfilter | lambda=400000")

df7$date = as.Date(df7$date)
ggplot(df7, aes(date, HP_trend, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price index 2010 = 100",
       subtitle = "Trend - HPfilter | lambda=400000")
