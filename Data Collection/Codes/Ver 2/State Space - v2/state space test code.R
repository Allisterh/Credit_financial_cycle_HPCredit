library(tidyverse)
library(tictoc)
library(ucminf)
library(numDeriv)

# Replicates Table 3 from Morley, 2007 JMCB
#
#
rm(list = ls())
# setwd("put working directory here")
setwd("D:/GitHub/HPCredit/Data Collection/Codes/Ver 2/State Space - v2") 

data_im <- read.table("D:/GitHub/HPCredit/Data Collection/MergedData-Raw.txt", header=TRUE, sep=",")
data_im = na.omit(data_im)

data_im <- data_im %>%
  filter(ID=="US")

data <- cbind(data_im$HPIndex,data_im$HHCredit)

source("trans.R") # Parameter constraints
source("lik_fcn.R") # Negative log likelihood function
source("filter_fcn.R") # Filter function

data = na.omit(data)
y <- 100*log(data)

T <- nrow(y)
#T <-49

START <- 2

prior <- 100

#=========================================================================#
# Maximum Likelihood Estimation
#=========================================================================#

# Initial values for optimisation routine

prmtr_in = c(-3.08,-24.92,-0.011,0.24649,
             -0.2259,0.53124,1.98897,-3.2913,
             0.0028,0.0005,
             -1.236,0.8003,9.4251,4.6742,
             -1.052,0.84236)

prmtr_in = runif(16, min=-1, max=1)

prmtr_in = t(prmtr_in)

trans(prmtr_in)

tic("ucminf")
# Initial paramter values
model = ucminf(prmtr_in,lik_fcn,hessian = TRUE,control = list(maxeval = 3000))
# Returns paramter estimates, -LL value, code
toc()

model$hessian
solve(model$hessian)
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
results = file("results.txt")

# Final Output
writeLines(c("Likelihood value is ", -model$value, 
             "code ", model$convergence, "",
             "Estimated parameters are:", c(t(prm_fnl),t(sd_fnl)), "",
             "Pre-transformed estimates are:", model$par,"",
             "Starting values:", prmtr_in), results)
close(results)





#=========================================================================#
# Impulse Response Functions
#=========================================================================#

filter_out = filter_fcn(model$par)
data = filter_out[[1]]
forcst = filter_out[[2]]

# Creates output file to store filtered dataset
write.csv(cbind(data[,1],data[,3],data[,5],forcst[,1:2]),"uc_yc.txt")

phi_y1 = prm_fnl[1]
phi_y2 = prm_fnl[2]
phi_c1 = prm_fnl[3]
phi_c2 = prm_fnl[4]

f_y = matrix(c(phi_y1,1,phi_y2,0),2,2)
f_c = matrix(c(phi_c1,1,phi_c2,0),2,2)

irf_fnl = c()
irf = 1
psi_ll = 0
psi_l = 1

for(j in 1:40){
  psi_t = phi_y1*psi_l + phi_y2*psi_ll
  irf = rbind(irf,psi_t)
  psi_ll = psi_l
  psi_l = psi_t
}

irf_fnl = cbind(irf_fnl,irf)
irf = 1
psi_ll = 0
psi_l = 1

for(j in 1:40){
  psi_t = phi_c1*psi_l + phi_c2*psi_ll
  irf = rbind(irf,psi_t)
  psi_ll = psi_l
  psi_l = psi_t
}

irf_fnl = cbind(irf_fnl,irf)

hlp = 0.5*matrix(1,nrow(irf_fnl),1) # Half Lives
hlm = -0.5*matrix(1,nrow(irf_fnl),1)

# Creates output file to store eigenvalues
eig_y = eigen(f_y)
eig_c = eigen(f_c)
eigenvalues = file("model1.eig.txt")
writeLines(c("Eigenvalues:",eig_y$values, abs(eig_y$values),
             eig_c$values, abs(eig_c$values)), eigenvalues)
close(eigenvalues)

# Creates output file to store irf dataset
write.csv(irf_fnl,"uc_yc_irf.txt");

#=========================================================================#
# Figures
#=========================================================================#
  
  par(mfrow=c(3,2))
  data_vec = seq(from = 1,to = T)
  
  # Income plot
  plot(data_vec[START:T],y[START:T,1],type = "n",main = "Income",
       xlab = "Month",ylab = "")
  lines(data_vec[START:T],y[START:T,1])
  
  # Consumption plot
  plot(data_vec[START:T],y[START:T,2],type = "l",main = "Consumption",
       xlab = "Month",ylab = "")
  
  # Permanent income plot
  plot(data_vec,data[,5],type = "l",main = "Permanent Income",
       xlab = "Month",ylab = "")
  
  # Permanent consumption plot
  plot(data_vec,data[,5] + prm_fnl[12],type = "l",main = "Permanent Consumption",
       xlab = "Month",ylab = "")
  
  # Transitory income plot
  plot(data_vec,data[,1],type = "l", main = "Transitory Income",
       xlab = "Month",ylab = "")
  lines(data_vec,matrix(0,T,1),lty = 2)
  
  # Transitory consumption plot
  plot(data_vec,data[,3],type = "l", main = "Transitory Consumption",
       xlab = "Month",ylab = "")
  lines(data_vec,matrix(0,T,1),lty = 2)
  
  
  par(mfrow=c(2,1))
  irf_vec = seq(from = 1,to = nrow(irf_fnl))

  # Plot income impulse response functions
  plot(irf_vec,irf_fnl[,1],type = "l", main = "IRF",
       xlab = "Quarter",ylab = "",lty = 1)
  lines(irf_vec,hlp,lty = 2)
  lines(irf_vec,matrix(0,nrow(irf_fnl)),lty = 3)
  lines(irf_vec,hlm,lty = 4)
  
  # Plot consumption impulse response functions
  plot(irf_vec,irf_fnl[,2],type = "l", main = "IRF",
       xlab = "Quarter",ylab = "",lty = 1)
  lines(irf_vec,hlp,lty = 2)
  lines(irf_vec,matrix(0,nrow(irf_fnl)),lty = 3)
  lines(irf_vec,hlm,lty = 4)
  