library(reshape)

#Dummy Creation
setwd("~/GitHub/HPCredit/Data Collection")
df <- data.frame(FC = rep(0,85),
                 FC_dummy= rep(0,85))
                 

seq0 = seq(as.Date("1999/1/1"), as.Date("2020/1/1"), "quarter")
df$date = seq0

seq0

write.table(df, "FCdummy1.csv", sep=',', row.names = FALSE)


####Import from excel
FC=read.csv("FCdummy.csv", sep=',')
FC1 = melt(FC, id.vars = "date", measure.vars = names(FC)[-1])
names(FC1)[2] = "ID"
names(FC1)[3] = "FCdummy"

write.table(FC1, "FCdummy.txt", sep=',', row.names = FALSE)
