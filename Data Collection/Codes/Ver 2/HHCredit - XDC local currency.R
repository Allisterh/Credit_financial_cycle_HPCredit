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
myvars <- c("borrowers_cty", "borrowers_country", "date", "obs_value")
df <- rates_plot[myvars]

names(df)[1]<-"ID"
names(df)[4]<-"HHCreditXDC"

write.table(df, "HHCreditXDC.txt", sep=",")

#-------------------------
#2. Data manipulation

df$ID = as.factor(df$ID)

#Trimming and cleaning data can be done here


#-------------------------
#3. Data inference
#Extract trends and cycles series from data

names(df)[4] = "value"

df7 <- df %>% group_by(ID) %>% 
  pdata.frame(., index = c("ID","date")) %>%
  mutate(HHCredit_localXDC_trend_1600 = mFilter::hpfilter(value, type = "lambda", freq = 1600)$trend)%>%
  mutate(HHCredit_localXDC_cycle_1600 = value - HHCredit_localXDC_trend_1600) %>%
  mutate(HHCredit_localXDC_trend_400k = mFilter::hpfilter(value, type = "lambda", freq = 400000)$trend)%>%
  mutate(HHCredit_localXDC_cycle_400k = value - HHCredit_localXDC_trend_400k)


names(df)[4] = "HHcredit"
df7$date = as.Date(df7$date)

write.table(df7, "HHCredit_localXDC_HPfilter.txt", sep=',' )

#names(df7)[4] = "HHCredit_localXDC_trend"
#names(df7)[5] = "HHCredit_localXDC_cycle"
#names(df7)[7] = "HHCredit_localXDC_gap"

#Graph the extracted series

#Cycle
ggplot(df7, aes(date, HHCredit_localXDC_cycle_1600, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Credit as local currency",
       subtitle = "cycle - HP filter decomp | lambda=1600")
ggsave("./graphs/hhcredit_localXDC_cycle.pdf", width=11, height=8.5)

#Trend
ggplot(melt(df7[,c(1:5)], c(1,2,3)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Credit as local currency",
       subtitle = "Trend - HP filter decomp | lambda=1600")
ggsave("./graphs/hhcredit_localXDC_trend.pdf", width=11, height=8.5)


#This part of code is to shape series into one graphs
df7$date = as.Date(df7$date)


df8 <- df7[,c(1,2,3,4,5,7)]
df9 <- df7[,c(1,2,3,6,8)]

df8 <- melt(df8, c(1,2,3))
df9 <- melt(df9, c(1,2,3))

#head(df8)
#table(df8$variable)

#Trend
ggplot(df8, aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Credit as local currency",
       subtitle = "Trend  vs observed value | lambda=1600 | 400k")
ggsave("./graphs/hhcredit_localXDC_trend_compare.pdf", width=11, height=8.5)


#Cycle
ggplot(df9, aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~borrowers_country) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "Household Credit as local currency",
       subtitle = "Cycle component | lambda=1600 | 400k")


ggsave("./graphs/hhcredit_localXDC_cycle_compare.pdf", width=11, height=8.5)
