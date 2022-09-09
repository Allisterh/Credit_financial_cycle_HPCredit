rm(list=ls())
#----
#1. Parameter Setup
country = "JP"


# Libraries
library(dplyr)
library(ggplot2)
library(zoo)
library(BIS)
library(reshape2)
library(rstudioapi)


setwd(dirname(getActiveDocumentContext()$path))
setwd("../../1.Latest/Paper2")
#---------------------
#1. Data Collection
#1.a. Credit 
#-------------------------



#Set up definition for get dataset function
datasets <- BIS::get_datasets()

#All avaiable country list (with HP available) 
#clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|SA|RS|SK|SI|SE|TR|US|ZA)"

#Country 

clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|RS|SK|SI|SE|TR|US)"
rates <- get_bis(datasets$url[datasets$name == "Credit-to-GDP gaps"], quiet = TRUE)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl(clist, borrowers_cty))%>%
  filter(grepl("^(A)", cg_dtype))


#filter(grepl("^(C|A)", cg_dtype)) 
# for both BIS generated Cycle and credit data

# 770 : Percentage of GDP
# USD and XDC is also available
# 
# Graph raw data for each country
# ggplot(rates_plot, aes(date, obs_value, color = borrowers_country)) +
#   geom_hline(yintercept = 0, linetype = "dashed",
#              color = "grey70", size = 0.02) +
#   geom_line(show.legend = FALSE) +
#   facet_wrap(~borrowers_country) +
#   theme_light() +
#   theme(panel.grid = element_blank()) +
#   labs(x = NULL, y = NULL,
#        title = "Total-Credit-to-non-financial",
#        subtitle = "as percentage of GDP")


#Saving Data
myvars <- c("borrowers_cty", "borrowers_country", "cg_dtype", "credit_gap_data_type", "date", "obs_value")
df <- rates_plot[myvars]

names(df)[1]<-"ID"

# df <-subset(df, date>as.Date(startdate))
# df <-subset(df, date<as.Date(enddate))

# write.table(df, "HHCredit.txt", sep=",")

#US code
df_1 = df %>%
  filter(ID == country)
df_1$date = as.Date(df_1$date)
filepath = sprintf("credit_%s.txt",country)
write.table(df_1, filepath, sep=',' )


plot(df_1$obs_value)
summary(df_1$obs_value)
summary(df_1$date)
