library(dplyr)
library(ggplot2)
library(zoo)
library(BIS)

setwd("~/GitHub/HPCredit/Data Collection")

datasets <- BIS::get_datasets()
#head(datasets, 20)

#country list
clist = "^(AT|AU|BE|BR|CA|CH|CL|CN|CO|CR|CZ|DE|DK|ES|EE|FI|FR|GB|GR|HU|ID|IN|IE|IS|IL|IT|JP|KR|LT|LU|LV|MX|NL|NO|NZ|PL|PT|RO|RU|SA|RS|SK|SI|SE|TR|US|ZA)"

rates <- get_bis(datasets$url[datasets$name == "Credit to the non-financial sector"], quiet = TRUE)

table(rates$tc_borrowers)
table(rates$borrowing_sector)

rates_plot <- rates %>%
  mutate(date = as.Date(as.yearqtr(date, "%Y-q%q"))) %>%
  filter(grepl(clist, borrowers_cty))%>%
  filter(grepl("^(C)", tc_borrowers))%>% #Private Credit
  filter(grepl("^(All sectors)", lending_sector))%>%
  filter(grepl("^(XDC)", unit_type))%>%
  #filter(grepl("^(A)", tc_adjust))%>%
  group_by(borrowers_cty) %>%
  mutate(growth = c(NA,diff(obs_value, lag=4))/obs_value)

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
       title = "Private Credit",
       subtitle = "Year to Year % change")


#Saving Data

myvars <- c("borrowers_cty", "date", "growth")
df <- rates_plot[myvars]

names(df)[1]<-"ID"
names(df)[3]<-"PrCredit"

df2 = df[df[,"date"]>="1999-01-01",]

write.table(df2, "PrCredit.txt", sep=",", row.names = FALSE)
