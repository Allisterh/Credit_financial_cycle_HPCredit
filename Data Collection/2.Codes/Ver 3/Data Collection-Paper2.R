rm(list=ls())
#----
# Set parameter
country = "AU"
# output = "credit_country.txt"
#----

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

# clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|RS|SK|SI|SE|TR|US)"
rates <- get_bis(datasets$url[datasets$name == "Credit to the non-financial sector"], quiet = TRUE)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(borrowers_cty == country)%>%
  filter(grepl("^(P)", tc_borrowers))%>%
  filter(grepl("^(A)", tc_lenders))%>%
  filter(grepl("^(A)", tc_adjust))%>%
  filter(grepl("^(770)", unit_type))

#770 : Percentage of GDP
#USD and XDC is also available
# 
#Graph raw data for each country
# ggplot(rates_plot, aes(date, obs_value, color = borrowers_country)) +
#   geom_hline(yintercept = 0, linetype = "dashed",
#              color = "grey70", size = 0.02) +
#   geom_line(show.legend = FALSE) +
#   facet_wrap(~borrowers_country) +
#   theme_light() +
#   theme(panel.grid = element_blank()) +
#   labs(x = NULL, y = NULL,
#        title = "Credit to household",
#        subtitle = "as percentage of GDP")


#Saving Data
myvars <- c("borrowers_cty", "borrowers_country", "date", "obs_value")
df <- rates_plot[myvars]

names(df)[1]<-"ID"
names(df)[4]<-"credit"

# df <-subset(df, date>as.Date(startdate))
# df <-subset(df, date<as.Date(enddate))

# write.table(df, "HHCredit.txt", sep=",")

#Export output data file
df_1 = df %>%
  filter(ID == country)
df_1$date = as.Date(df_1$date)
filepath = sprintf("credit_%s.txt",country)
write.table(df_1, filepath, sep=',' )


#-----------------------
#1.b Housing Prices
#-----------------------

#Preparing Dataset
detach("package:OECD", unload = TRUE)
get_datasets = BIS::get_datasets
datasets <- get_datasets()
head(datasets, 20)

#country list
#clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|SA|RS|SK|SI|SE|TR|US|ZA)"

#modified clist to match HP from BIS
# clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|FI|FR|GB|GR|HU|IE|ID|IN|IL|IT|JP|KR|LU|MX|NL|NO|NZ|PL|PT|RU|SA|SE|TR|US)"

rates <- get_bis(datasets$url[datasets$name == "Property prices: selected series"], quiet = TRUE)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(ref_area == country)%>%
  filter(grepl("^(R)", value))%>% #real values (vs nominal N)
  filter(grepl("^(628)", unit_measure)) 

df <- rates_plot[,c("ref_area", "reference_area", "date", "obs_value")]
names(df)[1]="ID"
names(df)[4]="HPIndex"

#771 vs 628
#628 : 2010 index = 100
#771 : Percentage change year on year

#table(rates_plot$ref_area)
#table(rates_plot$date)
#table(rates_plot$freq)
#table(rates$value)
#table(rates_plot$unit_of_measure)


#US code
df_1 = df %>%
  filter(ID == country)
df_1$date = as.Date(df_1$date)
filepath = sprintf("HPI_%s.txt",country)
write.table(df_1, filepath, sep=',' )

