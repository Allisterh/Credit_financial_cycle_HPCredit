library(dplyr)
library(ggplot2)
library(zoo)
library(BIS)

setwd("~/GitHub/HPCredit/Data Collection")

#Preparing Dataset
detach("package:OECD", unload = TRUE)
get_datasets = BIS::get_datasets
datasets <- get_datasets()
head(datasets, 20)

#country list
clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|SA|RS|SK|SI|SE|TR|US|ZA)"

#-----------------------
#House Price
#-----------------------

rates <- get_bis(datasets$url[datasets$name == "Property prices: selected series"], quiet = TRUE)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl(clist, ref_area))%>%
  filter(grepl("^(R)", value))%>%
  filter(grepl("^(Index, 2010 = 100)", unit_of_measure)) #771 vs 628

#table(rates_plot$ref_area)
#table(rates_plot$date)
#table(rates_plot$freq)
#table(rates_plot$value)
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


