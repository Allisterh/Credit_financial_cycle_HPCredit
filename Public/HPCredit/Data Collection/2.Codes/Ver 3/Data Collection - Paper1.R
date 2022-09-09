rm(list = ls())
library(dplyr)
library(ggplot2)
library(zoo)
library(reshape2)
library(rio)
library(data.table)
library(stringr)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
setwd("../../1.Latest/Paper1")

#-------
#Resources:
#https://www.bis.org/statistics/full_tc_csv.zip
#https://www.bis.org/statistics/full_spp_csv.zip

#---------------------
#1. Data Collection
#1.a. Household Credit
#-------------------------

# create a temporary directory
td <- tempdir()
# create a temporary file
tf <- tempfile(tmpdir = td, fileext = ".zip")
# download file from internet into temporary location
download.file("https://www.bis.org/statistics/full_tc_csv.zip", tf)
# list zip archive
file_names <- unzip(tf, list = TRUE)
# extract files from zip file
unzip(tf, exdir = td, overwrite = TRUE)
# use when zip file has only one file
data <- import(file.path(td, file_names$Name[1]))
# use when zip file has multiple files
#data_multiple <- lapply(file_names$Name, function(x) import(file.path(td, x)))
# delete the files and directories
rm(td)
rm(tf)


## Process data
rates_plot <- data %>%
  filter(.data[["Borrowing sector"]] == "H:Households & NPISHs") %>%
  filter(.data[["Unit type"]]
    == "770:Percentage of GDP")
  #filter(grepl(clist2, BORROWERS_CTY))%>%

#table(rates_plot[["Valuation method"]])


t_rates_plot <- transpose(rates_plot)

# get row and colnames in order
colnames(t_rates_plot) <- rownames(rates_plot)
rownames(t_rates_plot) <- colnames(rates_plot)
colnames(t_rates_plot) <- t_rates_plot[2, ] # changes col names to country id
#t_rates_plot <- t_rates_plot[-c(1:2,4:11),]
t_rates_plot <- t_rates_plot[-c(1:11), ] # Removing other naming scheme rows
t_rates_plot$date <- rownames(t_rates_plot)
#extract rownames -> column for date

dim(t_rates_plot)

### data type transforming
df <- t_rates_plot

df$date <- as.Date(as.yearqtr(df$date,           # Convert dates to quarterly
  format = "%Y-Q%q"))
#df <- xts(df, order.by=df$date)
#df <- ts_ts(df)

df1 <- df
latest_date <- file_names$Date[1]
name11 <- str_sub(names(df1[-ncol(df1)]), end = 2)
name21 <- str_sub(names(df1[-ncol(df1)]), start = 4)

colnames(df1) <- c(name11, "date")
colnames(df1)[which(names(df1) == "GB")] <- "UK"






#---------------------
#1. Data Collection
#1.b. House Price Index
#-------------------------

# create a temporary directory
td <- tempdir()
# create a temporary file
tf <- tempfile(tmpdir = td, fileext = ".zip")
# download file from internet into temporary location
download.file("https://www.bis.org/statistics/full_spp_csv.zip", tf)
# list zip archive
file_names <- unzip(tf, list = TRUE)
# extract files from zip file
unzip(tf, exdir = td, overwrite = TRUE)
# use when zip file has only one file
data <- import(file.path(td, file_names$Name[1]))
# use when zip file has multiple files
#data_multiple <- lapply(file_names$Name, function(x) import(file.path(td, x)))
# delete the files and directories
rm(td)
rm(tf)


## Process data
rates_plot <- data %>%
  filter(.data[["Unit of measure"]] == "Index, 2010 = 100") %>%
  filter(.data[["Value"]]
    == "Real")
  #filter(grepl(clist2, BORROWERS_CTY))%>%

#table(rates_plot[["Valuation method"]])


t_rates_plot <- transpose(rates_plot)

# get row and colnames in order
colnames(t_rates_plot) <- rates_plot$REF_AREA
rownames(t_rates_plot) <- colnames(rates_plot)
colnames(t_rates_plot) <- t_rates_plot[3, ] # changes col names to country id
#t_rates_plot <- t_rates_plot[-c(1:2,4:11),]
t_rates_plot <- t_rates_plot[-c(1:12), ] # Removing other naming scheme rows
t_rates_plot$date <- rownames(t_rates_plot)
#extract rownames -> column for date

dim(t_rates_plot)

### data type transforming
df <- t_rates_plot

df$date <- as.Date(as.yearqtr(df$date,           # Convert dates to quarterly
  format = "%Y-Q%q"))
#df <- xts(df, order.by=df$date)
#df <- ts_ts(df)

df2 <- df
latest_date <- file_names$Date[1]
#name12 <- str_sub(names(df1[-ncol(df1)]), end = 2)
#name22 <- str_sub(names(df1[-ncol(df1)]), start = 4)

#colnames(df2) <- c(name12, "date")
colnames(df2)[which(names(df2) == "GB")] <- "UK"

df11 <- df1
df22 <- df2


#------------
# 1.c Merge Data
name3 <- sort(Reduce(intersect, list(names(df1), names(df2))))
name3 <- name3[-which(name3 == "date")]
name4 <- c("4T", "5R", "XM")

for (i in seq_len(length(name4))) {
  name3 <- name3[-which(name3 == name4[i])]
} ## Remove "4T", "5R", "XM"

df1 <- df1[, c(name3, "date")]
df1$date <- as.Date(df1$date)
df1 <- reshape2::melt(df1, id = "date")
df2 <- df2[, c(name3, "date")]
df2$date <- as.Date(df2$date)
df2 <- reshape2::melt(df2, id = "date")

names(df1) <- c("date", "ID", "Credit")
names(df2) <- c("date", "ID", "HPI")

df1 <- na.omit(df1)
df2 <- na.omit(df2)

df <- merge(df1, df2, by = c("ID", "date"), all = TRUE)
df <- na.omit(df)


###----------
# Export data
filepath <- "fullsampleCreditHPI.csv"
write.table(df, filepath, sep = ",", row.names = FALSE)

sumdf <- df %>%
  group_by(ID) %>%
  summarize(min = min(date)) %>%
  arrange(min)

filepath <- "listofCountries.csv"
write.table(sumdf, filepath, sep = ",", row.names = FALSE)

sumdf <- sumdf %>%
  filter(min <= as.Date("1990-10-01"))
name3 <- as.character(sumdf[, "ID"][[1]])
filepath <- "shortlistofCountries.csv"
write.table(name3, filepath, sep = ",", row.names = FALSE)

## Export individual country
for (i in seq_len(length(name3))) {
  df3 <- df %>%
    filter(ID == name3[i]) %>%
    filter(date >= as.Date("1990-10-01"))
  filepath <- sprintf("creditHPI_%s.csv", name3[i])
  write.table(df3, filepath, sep = ",", row.names = FALSE)
}