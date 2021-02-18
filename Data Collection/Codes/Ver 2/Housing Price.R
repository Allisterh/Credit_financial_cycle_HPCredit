rm(list=ls())
library(dplyr)
library(ggplot2)
library(zoo)
library(BIS)
library(reshape2)

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
df_US = df %>%
  filter(ID == "US")

df_US$HPIndex_log = 100*log(df_US$HPIndex)
df_US$HPIndex_HPtrend = mFilter::hpfilter(df_US$HPIndex_log, type = "lambda", freq = 1600)$trend
df_US$HPIndex_HPcycle = mFilter::hpfilter(df_US$HPIndex_log, type = "lambda", freq = 1600)$cycle

df_US$date = as.Date(df_US$date)
write.table(df_US, "HPindex_HPfilter_US.txt", sep=',' )

#GB code
df_GB = df %>%
  filter(ID == "GB")

df_GB$HPIndex_log = 100*log(df_GB$HPIndex)
df_GB$HPIndex_HPtrend = mFilter::hpfilter(df_GB$HPIndex_log, type = "lambda", freq = 1600)$trend
df_GB$HPIndex_HPcycle = mFilter::hpfilter(df_GB$HPIndex_log, type = "lambda", freq = 1600)$cycle

df_GB$date = as.Date(df_GB$date)
write.table(df_GB, "HPindex_HPfilter_GB.txt", sep=',' )

#DE code
df_DE = df %>%
  filter(ID == "DE")

df_DE$HPIndex_log = 100*log(df_DE$HPIndex)
df_DE$HPIndex_HPtrend = mFilter::hpfilter(df_DE$HPIndex_log, type = "lambda", freq = 1600)$trend
df_DE$HPIndex_HPcycle = mFilter::hpfilter(df_DE$HPIndex_log, type = "lambda", freq = 1600)$cycle

df_DE$date = as.Date(df_DE$date)
write.table(df_DE, "HPindex_HPfilter_DE.txt", sep=',' )

#FR code
df_FR = df %>%
  filter(ID == "FR")

df_FR$HPIndex_log = 100*log(df_FR$HPIndex)
df_FR$HPIndex_HPtrend = mFilter::hpfilter(df_FR$HPIndex_log, type = "lambda", freq = 1600)$trend
df_FR$HPIndex_HPcycle = mFilter::hpfilter(df_FR$HPIndex_log, type = "lambda", freq = 1600)$cycle

df_FR$date = as.Date(df_FR$date)
write.table(df_FR, "HPindex_HPfilter_FR.txt", sep=',' )

#JP code
df_JP = df %>%
  filter(ID == "JP")

df_JP$HPIndex_log = 100*log(df_JP$HPIndex)
df_JP$HPIndex_HPtrend = mFilter::hpfilter(df_JP$HPIndex_log, type = "lambda", freq = 1600)$trend
df_JP$HPIndex_HPcycle = mFilter::hpfilter(df_JP$HPIndex_log, type = "lambda", freq = 1600)$cycle

df_JP$date = as.Date(df_JP$date)
write.table(df_JP, "HPindex_HPfilter_JP.txt", sep=',' )

#KR code
df_KR = df %>%
  filter(ID == "KR")

df_KR$HPIndex_log = 100*log(df_KR$HPIndex)
df_KR$HPIndex_HPtrend = mFilter::hpfilter(df_KR$HPIndex_log, type = "lambda", freq = 1600)$trend
df_KR$HPIndex_HPcycle = mFilter::hpfilter(df_KR$HPIndex_log, type = "lambda", freq = 1600)$cycle

df_KR$date = as.Date(df_KR$date)
write.table(df_KR, "HPindex_HPfilter_KR.txt", sep=',' )
#---------------------------------
#GRAPHING

#Cycle
ggplot(df_US, aes(date, HPIndex_HPcycle, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price Index, 2010 = 100",
       subtitle = "cycle - HP filter decomp | lambda=1600")
ggsave("./graphs/HPIndex_cycle.pdf", width=11, height=8.5)


#Trend
ggplot(df_US, aes(date, HPIndex_HPtrend, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price Index, 2010 = 100",
       subtitle = "Trend - HP filter decomp | lambda=1600")
ggsave("./graphs/HPIndex_cycle.pdf", width=11, height=8.5)


#Trend
ggplot(melt(df_US[,c(1,2,3,5,7)], c(1,2,3)), aes(date, HPIndex, color = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price Index, 100*log(value)",
       subtitle = "Trend - HP filter decomp | lambda=1600")
ggsave("./graphs/HPIndex_trend.pdf", width=11, height=8.5)


#cycle
ggplot(melt(df7[,c(1:5)], c(1,2,3)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price Index, 100*log(value)",
       subtitle = "Trend - HP filter decomp | lambda=1600")
ggsave("./graphs/HPIndex_trend.pdf", width=11, height=8.5)

#This part of code is to shape series into one graphs, show lambda=400k


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
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price Index, 2010 = 100",
       subtitle = "Trend  vs observed value | lambda=1600 | 400k")
ggsave("./graphs/HPIndex_trend_compare.pdf", width=11, height=8.5)


#Cycle
ggplot(df9, aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price Index, 2010 = 100",
       subtitle = "Cycle component | lambda=1600 | 400k")


ggsave("./graphs/HPIndex_cycle_compare.pdf", width=11, height=8.5)

