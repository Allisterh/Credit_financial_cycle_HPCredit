##
# Import 3 regression tables
  #list out links to the 3 files/ output 
rm(list=ls())

country = 'GB'

working_dir="D:/GitHub/HPCredit/Regression"
setwd(working_dir) 

filepath = sprintf("VAR_2/Output/Reg_%s.csv",country)
df1 <- read.table(filepath, header=TRUE, sep=",")

filepath = sprintf("VAR_2_crosscycle/Output/Reg_%s.csv",country)
df2 <- read.table(filepath, header=TRUE, sep=",")

filepath = sprintf("VAR_2_crosscycle_1stlagonly/Output/Reg_%s.csv",country)
df3 <- read.table(filepath, header=TRUE, sep=",")


# Merge 3 df

# Clean up table

# Export table file