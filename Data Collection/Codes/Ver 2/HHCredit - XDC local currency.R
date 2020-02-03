library(dplyr)
library(ggplot2)
library(zoo)
library(BIS)

setwd("D:/GitHub/HPCredit/Data Collection/")

#---------------------
#1. Data Collection

#Set up definition for get dataset function
datasets <- BIS::get_datasets()

#All avaiable country list (with HP available) 
#clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|SA|RS|SK|SI|SE|TR|US|ZA)"

#Country 

clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|RS|SK|SI|SE|TR|US)"
rates <- get_bis(datasets$url[datasets$name == "Credit to the non-financial sector"], quiet = TRUE)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl(clist, borrowers_cty))%>%
  filter(grepl("^(H)", tc_borrowers))%>%
  filter(grepl("^(All sectors)", lending_sector))%>%
  filter(grepl("^(XDC)", unit_type))

#770 : Percentage of GDP
#USD and XDC is also available

#Graph raw data for each country
ggplot(rates_plot, aes(date, obs_value, color = borrowers_country)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Credit to household",
       subtitle = "local currency")


#Saving Data
myvars <- c("borrowers_cty", "date", "obs_value")
df <- rates_plot[myvars]

names(df)[1]<-"ID"
names(df)[3]<-"HHCredit"

write.table(df, "HHCredit.txt", sep=",")

#-------------------------
#2. Data manipulation

df$ID = as.factor(df$ID)

#Trimming and cleaning data can be done here


#-------------------------
#3. Data inference
#Extract trends and cycles series from data
library(plm)

names(df)[1] = "ID"
names(df)[3] = "value"

df7 <- df %>% group_by(variable) %>% 
  pdata.frame(., index = c("ID","date")) %>%
  mutate(HHCredit_GDP_trend = mFilter::hpfilter(value, type = "lambda", freq = 1600)$trend)%>%
  mutate(HHCredit_GDP_cycle = mFilter::hpfilter(value, type = "lambda", freq = 1600)$cycle)%>%
  mutate(HHCredit_GDP_gap = value - mFilter::hpfilter(value, type = "lambda", freq = 1600)$trend)

write.table
write.table(df7, "HHCredit_GDP_HPfilter.txt", sep=',' )

names(df7)[1] = "ID"
names(df7)[2] = "date"

df7$date = as.Date(df7$date)
#names(df7)[4] = "HHCredit_GDP_trend"
#names(df7)[5] = "HHCredit_GDP_cycle"
#names(df7)[7] = "HHCredit_GDP_gap"

#Graph the extracted series

#Cycle
ggplot(df7, aes(date, HHCredit_GDP_cycle, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Credit - local currency",
       subtitle = "cycle - HP decomp | lambda=1600")

#Trend
ggplot(df7, aes(date, HHCredit_GDP_trend, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Credit - local currency",
       subtitle = "Trend - HP decomp | lambda=1600")

#Gap
ggplot(df7, aes(date, HHCredit_GDP_gap, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ID) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Credit - local currency",
       subtitle = "Gap = observed value - HP filter trend | lambda=1600")