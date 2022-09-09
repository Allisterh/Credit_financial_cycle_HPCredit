rm(list=ls())

#Import data
country = 'US'

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

VAR2_filepath = sprintf("Bayesian_UC_VAR2_drift/OutputData/Reg_%s.csv",country)
VAR2x1_filepath = sprintf("Bayesian_UC_VAR2_drift_Crosscycle1lag/OutputData/Reg_%s.csv",country)
VAR2x2_filepath = sprintf("Bayesian_UC_VAR2_drift_Crosscycle2lags/OutputData/Reg_%s.csv",country)

df1 <- read.table(VAR2_filepath, header=FALSE, sep=",")
df2 <- read.table(VAR2x1_filepath, header=FALSE, sep=",")
df3 <- read.table(VAR2x2_filepath, header=FALSE, sep=",")

df1 <- t(df1)
df2 <- t(df2)
df3 <- t(df3)
#Combine data
## SD table
df = data.frame(matrix(ncol = 6, nrow = 15))

df[1:2,1:2] = df1[1:2,1:2]
df[5:6,1:2] = df1[3:4,1:2]
df[9:15,1:2] = df1[5:11,1:2]


df[1:3,3:4] = df2[1:3,1:2]

df[7,3:4] = df2[4,1:2]
df[5:6,3:4] = df2[5:6,1:2]

df[9:15,3:4] = df2[7:13,1:2]


df[1:4,5:6] = df3[1:4,1:2]
df[5:6,5:6] = df3[7:8,1:2]
df[7:8,5:6] = df3[5:6,1:2]
df[9:15,5:6] = df3[9:15,1:2]

#Export csv
df4 <- read.table("symbols.csv", header=FALSE, sep=",")
df=cbind(df4,df)


filepath=sprintf('RegComb_%s.csv',country)
write.table(df, filepath,  
            na = "",
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            sep = " , ")
#            eol=" \\\\[2pt] \r\n",
#            quote = FALSE)


#install.packages("kableExtra")
library('kableExtra')
library(dplyr)
library(knitr)
rownames(df) <- df[,1]
df<-df[,-1]
colnames(df) <- c("Est", "SD", "Est", "SD", "Est", "SD")

options(knitr.kable.NA = '')

#df = df %>% mutate_if(is.numeric, format, digits=4)

#kbl(data.frame(x=rnorm(10), y=rnorm(10), x= rnorm(10)), digits = c(1, 4, 4))

kbl(df, digits = c(4, 4, 4, 4, 4, 4)) %>%
  kable_classic("striped") %>%
  add_header_above(c("Parameters" = 1, "VAR2" = 2, "VAR2 1-cross lag" = 2, "VAR2 2-cross lags" = 2)) %>%
  footnote(general="UK Bayesian regression results") %>%
  kable_styling(latex_options="scale_down")


# ## 10-90 percentile

# #Combine data
# ## percentile
# df = data.frame(matrix(ncol = 9, nrow = 15))

# df[1:2,1:3] = df1[1:2,c(1,3:4)]
# df[5:6,1:3] = df1[3:4,c(1,3:4)]
# df[9:15,1:3] = df1[5:11,c(1,3:4)]


# df[1:3,4:6] = df2[1:3,c(1,3:4)]

# df[7,4:6] = df2[4,c(1,3:4)]
# df[5:6,4:6] = df2[5:6,c(1,3:4)]

# df[9:15,4:6] = df2[7:13,c(1,3:4)]


# df[1:4,7:9] = df3[1:4,c(1,3:4)]
# df[5:6,7:9] = df3[7:8,c(1,3:4)]
# df[7:8,7:9] = df3[5:6,c(1,3:4)]
# df[9:15,7:9] = df3[9:15,c(1,3:4)]

# #Export csv
# df4 <- read.table("symbols.csv", header=FALSE, sep=",")
# df=cbind(df4,df)

# filepath=sprintf('RegCombPercentile_%s.csv',country)
# write.table(df, filepath,  
#             na = "",
#             row.names = FALSE,
#             col.names = FALSE,
#             append = FALSE,
#             sep = " , ")
# #            eol=" \\\\[2pt] \r\n",
# #            quote = FALSE)

# #install.packages("kableExtra")
# library('kableExtra')
# library(dplyr)
# library(knitr)
# rownames(df) <- df[,1]
# df<-df[,-1]
# colnames(df) <- c("Median", "10%", "90%", "Median", "10%", "90%", "Median", "10%", "90%")

# options(knitr.kable.NA = '')

# #df = df %>% mutate_if(is.numeric, format, digits=4)

# #kbl(data.frame(x=rnorm(10), y=rnorm(10), x= rnorm(10)), digits = c(1, 4, 4))

# kbl(df, digits = c(4, 4, 4, 4, 4, 4, 4, 4, 4)) %>%
#   kable_paper("striped") %>%
#   add_header_above(c("Parameters" = 1, "VAR2" = 3, "VAR2 1-cross lag" = 3, "VAR2 2-cross lags" = 3)) %>%
#   footnote(general="UK Bayesian regression results") %>%
#   kable_styling(latex_options="scale_down") %>%
#   column_spec(5:7, bold = c(0,0,1,0,0,0,1,0,0,0,0,0,0,0,0))




## 10-90 percentile

#Combine data
## percentile
df = data.frame(matrix(ncol = 9, nrow = 15))

df[1:2,1:3] = df1[1:2,c(1,3:4)]
df[5:6,1:3] = df1[3:4,c(1,3:4)]
df[9:15,1:3] = df1[5:11,c(1,3:4)]


df[1:3,4:6] = df2[1:3,c(1,3:4)]

df[7,4:6] = df2[4,c(1,3:4)]
df[5:6,4:6] = df2[5:6,c(1,3:4)]

df[9:15,4:6] = df2[7:13,c(1,3:4)]


df[1:4,7:9] = df3[1:4,c(1,3:4)]
df[5:6,7:9] = df3[7:8,c(1,3:4)]
df[7:8,7:9] = df3[5:6,c(1,3:4)]
df[9:15,7:9] = df3[9:15,c(1,3:4)]

library(tidyr)
library(dplyr)
df0 <- df
df = data.frame(matrix(ncol = 6, nrow = 15))

f <- function(x){ sprintf('%.4f',x)}
df0<- lapply(df0, f)

df0<-data.frame(df0)
df[,1] = df0[,1]
df[,2] <- df0[,2:3] %>% 
   unite(combined1, X2, X3 , sep = ", ", remove = TRUE) %>%
   mutate(combined1 = sprintf("[%s]", combined1))

df[,3] = df0[,4]
df[,4] <- df0[,5:6] %>% 
   unite(combined2, X5, X6 , sep = ", ", remove = TRUE) %>%
   mutate(combined2 = sprintf("[%s]", combined2))

df[,5] = df0[,7]
df[,6] <- df0[,8:9] %>% 
   unite(combined3, X8, X9 , sep = ", ", remove = TRUE) %>%
   mutate(combined3 = sprintf("[%s]", combined3))

#https://stackoverflow.com/questions/3357743/replacing-character-values-with-na-in-a-data-frame
df[ df == "[NA, NA]" ] <- NA
df[ df == "NA" ] <- NA

#Export csv
df4 <- read.table("symbols.csv", header=FALSE, sep=",")
df=cbind(df4,df)
df4 <- read.table("symbols_desc.csv", header=FALSE, sep=",")
df=cbind(df4,df)

filepath=sprintf('RegCombPercentile_%s.csv',country)
write.table(df, filepath,  
            na = "",
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            sep = " , ")
#            eol=" \\\\[2pt] \r\n",
#            quote = FALSE)

#install.packages("kableExtra")
library('kableExtra')
library(dplyr)
library(knitr)
#rownames(df) <- df[,1]
#df<-df[,-1]
colnames(df) <- c("Description","Parameters","Median", "[10%, 90%]", "Median", "[10%, 90%]", "Median", "[10%, 90%]")

options(knitr.kable.NA = '')

#df = df %>% mutate_if(is.numeric, format, digits=4)

#kbl(data.frame(x=rnorm(10), y=rnorm(10), x= rnorm(10)), digits = c(1, 4, 4))

kbl(df, digits = c(4, 4, 4, 4, 4, 4), row.names=FALSE) %>%
  kable_paper("striped") %>%
  add_header_above(c(" " = 2, "VAR2" = 2, "VAR2 1-cross lag" = 2, "VAR2 2-cross lags" = 2)) %>%
  footnote(general="UK Bayesian regression results") %>%
  kable_styling(latex_options="scale_down") %>%
  column_spec(5:7, bold = c(0,0,1,0,0,0,1,0,0,0,0,0,0,0,0))




## Cross-countries evidence

filepath <- "../Data Collection/1.Latest/Paper1/shortlistofCountries_full.csv"

name1 <- read.table(filepath, header = TRUE, sep = ",")

df = data.frame()

for (i in seq_len(nrow(name1))) {
country <- name1[i, 1]

filepath <- sprintf("../Regression/Bayesian_UC_VAR2_drift_Crosscycle1lag/OutputData/Reg_%s.csv", country)
df1 <- read.table(filepath, header = FALSE, sep = ",")
df1 <- df1[c(1,3:4),3:4]
df2 = data.frame(matrix(ncol = 7, nrow = 1))
df2[1,2:4] <- t(df1[1:3,1])
df2[1,5:7] <- t(df1[1:3,2])
df3 = data.frame(matrix(ncol = 5, nrow = 1))
df3[1,1] <- name1[i, 2]
df3[1,2] <- df2[1,2]
df3[1,4] <- df2[1,5]
f <- function(x){ sprintf('%.4f',x)}
df2<- data.frame(lapply(df2, f))
df2[1,1] <- name1[i, 2]
df3[1,3] <- df2[,3:4] %>% 
   unite(combined1, X3, X4 , sep = ", ", remove = TRUE) %>%
   mutate(combined1 = sprintf("[%s]", combined1))
df3[1,5] <- df2[,6:7] %>% 
   unite(combined1, X6, X7 , sep = ", ", remove = TRUE) %>%
   mutate(combined1 = sprintf("[%s]", combined1))
df <- rbind(df,df3)
}

colnames(df) <- c("Country","Median", "[10%, 90%]", "Median", "[10%, 90%]")

kbl(df, digits = c(4, 4, 4, 4, 4), row.names=FALSE) %>%
  kable_paper("striped") %>%
  add_header_above(c(" " = 1, "$\\phi^{x1}_y$" = 2, "$\\phi^{x1}_h$" = 2)) %>%
  footnote(general="Cross countries") %>%
  kable_styling(latex_options="scale_down")