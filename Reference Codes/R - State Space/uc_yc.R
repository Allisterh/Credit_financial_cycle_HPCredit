# Replicates Table 3 from Morley, 2007 JMCB
#
#
rm(list = ls())
# setwd("put working directory here")
library(tictoc)

setwd("D:/Github/HPCredit/Data Collection/Codes/Ver 2/State Space/")

library(numDeriv)
library(ucminf)

data_im <- read.table("yc.txt")
data <- cbind(data_im[9:214,1],data_im[9:214,3])

source("trans.R") # Parameter constraints
source("lik_fcn.R") # Negative log likelihood function
source("filter_fcn.R") # Filter function

y <- 100*log(data)

T <- nrow(y)

START <- 2

prior <- 100

#=========================================================================#
# Maximum Likelihood Estimation
#=========================================================================#

# Initial values for optimisation routine
prmtr_in = c(1,-4,0.74382,-5.07080, 
             0.51159,
             -0.25900, 0.41104,7.31927,
              1.02369,-1.50885,-0.76931,
             -0.16347,0.96781)
prmtr_in = t(prmtr_in)

trans(prmtr_in)


tic("ucminf model")
# Initial paramter values
model = ucminf(prmtr_in,lik_fcn,hessian = TRUE,control = list(maxeval = 3000))
# Returns paramter estimates, -LL value, code
toc()


# Final parameter values
prm_fnl = t(trans(model$par))

# Use Hessian to find parameter standard errors
hessn0 = model$hessian
cov0 = solve(hessn0)

grdn_fnl = jacobian(trans, model$par)
cov = grdn_fnl%*%cov0%*%t(grdn_fnl)
sd_fnl = sqrt(abs(diag(cov)))
sd_out = sqrt(abs(diag(cov0)))

# Create output file to store results
results = file("results.txt")
writeLines(c("Starting values:", prmtr_in),results)

# Final Output
writeLines(c("Likelihood value is ", -model$value, 
             "code ", model$convergence, "",
             "Estimated parameters are:", c(t(prm_fnl),t(sd_fnl)), "",
             "Pre-transformed estimates are:", model$par),results)
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
  