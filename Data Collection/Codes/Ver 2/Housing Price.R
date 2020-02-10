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
  filter(grepl("^(628)", unit_measure)) 

df <- rates_plot[,c("ref_area", "reference_area", "date", "obs_value")]
names(df)[1]="variable"
names(df)[4]="value"

#771 vs 628
#628 : 2010 index = 100
#771 : Percentage change year on year

#table(rates_plot$ref_area)
#table(rates_plot$date)
#table(rates_plot$freq)
#table(rates$value)
#table(rates_plot$unit_of_measure)


df7 <- df %>% group_by(variable) %>% 
  pdata.frame(., index = c("variable","date")) %>%
  mutate(HPIndex_trend_1600 = mFilter::hpfilter(value, type = "lambda", freq = 1600)$trend)%>%
  mutate(HPIndex_cycle_1600 = value - HPIndex_trend_1600) %>%
  mutate(HPIndex_trend_400k = mFilter::hpfilter(value, type = "lambda", freq = 400000)$trend)%>%
  mutate(HPIndex_cycle_400k = value - HPIndex_trend_400k)

names(df7)[1] = "ID"
names(df7)[4] = "HPIndex"
names(df7)[5] = "HPIndex_trend_1600"
names(df7)[6] = "HPIndex_cycle_1600"
names(df7)[7] = "HPIndex_trend_400k"
names(df7)[8] = "HPIndex_cycle_400k"

df7$date = as.Date(df7$date)

write.table(df7, "HPindex_HPfilter.txt", sep=',' )

#---------------------------------
#GRAPHING

#Cycle
ggplot(df7, aes(date, HPIndex_cycle_1600, color = ID)) +
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
ggplot(melt(df7[,c(1:5)], c(1,2,3)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  facet_wrap(~reference_area) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = "House Price Index, 2010 = 100",
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

