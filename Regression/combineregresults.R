 rm(list=ls())


#Import data
country = 'GB'
setwd("D:/GitHub/HPCredit/Paper/Regression")

VAR2_filepath = "D:/GitHub/HPCredit/Regression/VAR_2/Output/Reg_GB.csv"
VAR2x2_filepath = "D:/GitHub/HPCredit/Regression/VAR_2_crosscycle/Output/Reg_GB.csv"
VAR2x1_filepath = "D:/GitHub/HPCredit/Regression/VAR_2_crosscycle_1stlagonly/Output/Reg_GB.csv"

df1 <- read.table(VAR2_filepath, header=FALSE, sep=",")
df2 <- read.table(VAR2x1_filepath, header=FALSE, sep=",")
df3 <- read.table(VAR2x2_filepath, header=FALSE, sep=",")


#Combine data

df = data.frame(matrix(ncol = 6, nrow = 15))
df[1:13,5:6] = df3[1:13,]
df[15,5] = df3[14,1]
df[1:3,3:4] = df2[1:3,]
df[5:7,3:4] = df2[4:6,]
df[9:13,3:4] = df2[7:11,]
df[15,3] = df2[12,1]
df[1:2,1:2] = df1[1:2,]
df[5:6,1:2] = df1[3:4,]
df[9:15,1:2] = df1[5:11,]
df[15,2] = NA

#Export csv
df4 <- read.table("symbols.csv", header=FALSE, sep=",")
df=cbind(df4,df)

write.table(df, "RegComb_GB.csv",  
            na = "",
            row.names = FALSE,
            col.names = FALSE,
            append = FALSE,
            sep = " & ",
            eol=" \\\\[2pt] \r\n",
            quote = FALSE)

