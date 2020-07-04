#Load libraries
rm(list = ls())

library(tidyverse)
library(tictoc)
library(ucminf)
library(numDeriv)
#library(matrixcalc)
library(DataCombine)
library(dplyr)
library(reshape2)
library(ggplot2)

source("trans.R") # Parameter constraints
source("lik_fcn.R") # Negative log likelihood function
source("filter_fcn.R") # Filter function


# Replicates Table 3 from Morley, 2007 JMCB
#

# setwd("put working directory here")
setwd("D:/GitHub/HPCredit/Regression/Ver 2 - UC with random walk drift term/R") 

country = 'GB'
#automate file name  sprintf("input%s.txt",country)

#Read raw data
data_filepath = sprintf("../Data/Input/data_%s.txt",country)
data <- read.table(data_filepath, header=TRUE, sep=",")

#=========================================================================#
# UC model regression
#=========================================================================#

data = na.omit(data)
y <- 100*log(data)
START = 2

#Setting Priors

t_y_prior = y[1,1]
t_h_prior = y[1,2]
y=y[-1,]
T = nrow(y)

  #Loop for solvable hessian matrix

number = 1
det_model = TRUE
while (det_model==TRUE)
{
  #Setting priors
  sig_ty_prior = as.numeric(100+sample(1:100, 1))
  sig_gy_prior = as.numeric(sample(1:10, 1))
  sig_th_prior = as.numeric(100+sample(1:100, 1))
  sig_gh_prior = as.numeric(sample(1:10, 1))
  sig_tyth_prior = as.numeric(sample(1:10, 1))
  #randomize cross trend covariance to be large positive or negative number.
  # m = sample(1:2,1)-1
  # if (m == 0) {m = -1}
  # sig_tyth_prior = as.numeric(m*(50+sample(1:50, 1)))

  prior = c(t_y_prior, t_h_prior, sig_ty_prior, sig_gy_prior, sig_th_prior, sig_gh_prior, sig_tyth_prior)
  
  prmtr_in = runif(12, min=-20, max=20)
  
  # det_ft = prior_setting(prmtr_in,prior)
  # det_ft
  
  tic("ucminf")
  # Initial paramter values
  model = ucminf(prmtr_in,lik_fcn,hessian = TRUE,control = list(maxeval = 3000))
  # Returns paramter estimates, -LL value, code
  toc()
  
  print(paste("Number of loops: ", number))
  number = number + 1
  det_model=is.nan(det(model$hessian))
}

# Returns paramter estimates, -LL value, code

# Final parameter values
prm_fnl = t(trans(model$par))

model$par
prm_fnl

# Use Hessian to find parameter standard errors
hessn0 = model$hessian
cov0 = solve(hessn0)

grdn_fnl = jacobian(trans, model$par) #trans function used here
cov = grdn_fnl%*%cov0%*%t(grdn_fnl)
sd_fnl = sqrt(abs(diag(cov)))
sd_out = sqrt(abs(diag(cov0)))
 
# Create output file to store results
setwd("../Output") 
result_filepath = sprintf("R_result_%s.txt",country)
results = file(result_filepath)

# Final Output
writeLines(c("Likelihood value is ", -model$value, 
             "code ", model$convergence, "",
             "Estimated parameters are:", c(t(prm_fnl),t(sd_fnl)), "",
             "Pre-transformed estimates are:", model$par,"",
             "Starting values:", prmtr_in,"",
             "Staring priors:", prior,""), results)
close(results)
setwd("../R") 

#export priors and starting values
write.table(prmtr_in,sprintf("../Data/R_prmtr_in_%s.txt",country),col.names = FALSE, row.names = FALSE)
write.table(prior,sprintf("../Data/R_prior_%s.txt",country),col.names = FALSE, row.names = FALSE)

#write regression results to csv
reg = cbind(t(prm_fnl),matrix(sd_fnl,12,1))
reg = rbind(reg,c(-model$value,0))
write.table(reg,sprintf("../Output/Reg_%s.csv",country),col.names = FALSE, row.names = FALSE)


#=========================================================================#
# Forecast data
#=========================================================================#

filter_out = filter_fcn(model$par)
data = filter_out[[1]]
forcst = filter_out[[2]]

# Creates output file to store filtered dataset
write.table(cbind(data[,1],data[,3],data[,4],data[,6],forcst[,1:2]),sprintf("../Data/uc_yc_%s.txt",country),sep=',',row.names = FALSE,col.names = FALSE)


#=========================================================================#
# Impulse Response Functions
#=========================================================================#

# 
# df1 <- read.table( sprintf("../Output/Reg_%s.csv",country) , header=FALSE, sep=",")
# prm_fnl = df1[,1]

phi_y1 = prm_fnl[1]
phi_yh = prm_fnl[3]

phi_h1 = prm_fnl[4]
phi_hy = prm_fnl[6]


irf_fnl = c()
irf = 1
psi_ll = 0
psi_l = 1

for(j in 1:40){
  psi_t = phi_y1*psi_l + phi_yh*psi_l
  irf = rbind(irf,psi_t)
  psi_ll = psi_l
  psi_l = psi_t
}

irf_fnl = cbind(irf_fnl,irf)
irf = 1
psi_ll = 0
psi_l = 1

for(j in 1:40){
  psi_t = phi_h1*psi_l + phi_hy*psi_l
  irf = rbind(irf,psi_t)
  psi_ll = psi_l
  psi_l = psi_t
}

irf_fnl = cbind(irf_fnl,irf)

hlp = 0.5*matrix(1,nrow(irf_fnl),1) # Half Lives
hlm = -0.5*matrix(1,nrow(irf_fnl),1)

irf_vec = seq(from = 1,to = nrow(irf_fnl))

par(mar=c(1,1,1,1))


pdf(file= sprintf("../Output/Graphs/IRF_%s.pdf",country) )
par(mfrow=c(2,1))
# Plot income impulse response functions
plot(irf_vec,irf_fnl[,1],type = "l", main = "Credit IRF",
     xlab = "Quarter",ylab = "",lty = 1)
lines(irf_vec,hlp,lty = 2)
lines(irf_vec,matrix(0,nrow(irf_fnl)),lty = 3)
lines(irf_vec,hlm,lty = 4)

# Plot consumption impulse response functions
plot(irf_vec,irf_fnl[,2],type = "l", main = "Housing Price IRF",
     xlab = "Quarter",ylab = "",lty = 1)
lines(irf_vec,hlp,lty = 2)
lines(irf_vec,matrix(0,nrow(irf_fnl)),lty = 3)
lines(irf_vec,hlm,lty = 4)

dev.off()


#=========================================================================#
# Plotting UC and HP filter data
#=========================================================================#

HP_filepath = sprintf("../Data/Input/HPindex_HPfilter_%s.txt",country)
df1 <- read.table(HP_filepath, header=TRUE, sep=",")

# df1 <- na.omit(df1[-c(2)]) Remove country column

Credit_filepath = sprintf("../Data/Input/Credit_HPfilter_%s.txt",country)
df2 <- read.table(Credit_filepath, header=TRUE, sep=",")
df2 <- na.omit(df2[-c(2)]) #Remove country name column because redundancy

df <- merge(df1, df2, by=c("ID","date"))
df = df[-1,] #to account for prior quarter

df3 <- read.table(sprintf("../Data/uc_yc_%s.txt",country), header=FALSE, sep=",")

df = df[-nrow(df),]
df = cbind(df,df3)
df$date = as.Date(df$date)

#head(df8)
#table(df8$variable)


#Cycles var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "HPIndex_HPcycle", "V2", "V4")

#Credit Cycle var name list
varlist2 = c("ID", "date", "Credit_HPcycle", "V2")

df6 = df[varlist2]
names(df6)[4]="UC_Credit_Cycle"
names(df6)[3]="HP_Credit_Cycle"

ggplot(melt(df6, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Credit cycle: %s",country)  )
ggsave( sprintf("../Output/graphs/Credit_cycle_%s.pdf",country) , width=8, height=5)

#HP Cycle
varlist3 = c("ID", "date", "HPIndex_HPcycle", "V4")
df7 = df[varlist3]
names(df7)[4]="UC_HPI_Cycle"
names(df7)[3]="HP_HPI_Cycle"


ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Housing Price cycle: %s",country)  )
ggsave( sprintf("../Output/graphs/HP_cycle_%s.pdf",country) , width=8, height=5)

#Trends 
#Credit Trends

varlist4 = c("ID", "date", "Credit_HPtrend", "V1", "Credit_log")
df7 = df[varlist4]
names(df7)[4]="UC_Credit_Trend"
names(df7)[3]="HP_Credit_Trend"
names(df7)[5]="Series"
df7[c(3:5)] = exp(df7[c(3:5)]/100)

ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Credit Trend: %s , as percentage of GDP",country)  )
ggsave( sprintf("../Output/graphs/Credit_Trend_%s.pdf",country) , width=8, height=5)


#Housing Price Index Trends

varlist5 = c("ID", "date", "HPIndex_HPtrend", "V3", "HPIndex_log")
df7 = df[varlist5]
names(df7)[4]="UC_Credit_Trend"
names(df7)[3]="HP_Credit_Trend"
names(df7)[5]="Series"
df7[c(3:5)] = exp(df7[c(3:5)]/100)

ggplot(melt(df7, c(1,2)), aes(date, value, color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", size = 0.02) +
  geom_line(show.legend = TRUE) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL,
       title = sprintf("Housing Price Index Trend: %s , Index 2010=100",country)  )
ggsave( sprintf("../Output/graphs/HP_Trend_%s.pdf",country) , width=8, height=5)

