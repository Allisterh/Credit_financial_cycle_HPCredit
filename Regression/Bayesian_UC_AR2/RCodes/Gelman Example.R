rm(list=ls())

library(tidyverse)
library(mvtnorm)

# set seed for replicability
set.seed(8675309)

# create a N x k matrix of covariates
N = 250
K = 3

covariates = replicate(K, rnorm(n = N))
colnames(covariates) = c('X1', 'X2', 'X3')

# create the model matrix with intercept
X = cbind(Intercept = 1, covariates)

# create a normally distributed variable that is a function of the covariates
coefs = c(5, .2, -1.5, .9)
sigma = 2
mu = X %*% coefs
y  = rnorm(N, mu, sigma)

# same as
# y = 5 + .2*X1 - 1.5*X2 + .9*X3 + rnorm(N, mean = 0, sd = 2)

# Run lm for later comparison; but go ahead and examine now if desired
fit_lm = lm(y ~ ., data = data.frame(X[, -1]))
# summary(fit_lm)

##---- 2nd part
log_posterior <- function(X, y, th) {
  # Args
  # X: the model matrix
  # y: the target vector
  # th: theta, the current parameter estimates
  
  beta   = th[-length(th)]            # reg coefs to be estimated
  sigma  = th[length(th)]             # sigma to be estimated
  sigma2 = sigma^2
  mu = X %*% beta
  
  # priors are b0 ~ N(0, sd=10), sigma2 ~ invGamma(.001, .001)
  priorbvarinv = diag(1/100, 4) 
  prioralpha   = priorbeta = .001
  
  if (is.nan(sigma) | sigma<=0) {     # scale parameter must be positive, so post
    return(-Inf)                      # density is zero if it jumps below zero
  }
  # log posterior in this conjugate setting. conceptually it's (log) prior +
  # (log) likelihood. (See commented 'else' for alternative)
  # else {
  #   -.5*nrow(X)*log(sigma2) - (.5*(1/sigma2) * (crossprod(y-mu))) +
  #     -.5*ncol(X)*log(sigma2) - (.5*(1/sigma2) * (t(beta) %*% priorbvarinv %*% beta)) +
  #     -(prioralpha + 1)*log(sigma2) + log(sigma2) - priorbeta/sigma2
  # }
  else {
    ll = mvtnorm::dmvnorm(y, mean=mu, sigma=diag(sigma2, length(y)), log=T)
    priorb = mvtnorm::dmvnorm(beta, mean=rep(0, length(beta)), sigma=diag(100, length(beta)), log=T)
    priors2 = dgamma(1/sigma2, prioralpha, priorbeta, log=T)
    logposterior = ll + priorb + priors2
    logposterior
  }
}

##---- 3rd Part
gradient_theta <- function(X, y, th) {
  d = length(th)
  e = .0001
  diffs = numeric(d)
  
  for (k in 1:d) {
    th_hi = th
    th_lo = th
    th_hi[k] = th[k] + e
    th_lo[k] = th[k] - e
    diffs[k] = (log_posterior(X, y, th_hi) - log_posterior(X, y, th_lo)) / (2 * e)
  }
  
  diffs
}

##---- 4th part

hmc_iteration <- function(X, y, th, epsilon, L, M) {
  # Args
  # epsilon: the stepsize
  # L: the number of leapfrog steps
  # M: a diagonal mass matrix
  
  # initialization
  M_inv = 1/M
  d   = length(th)
  phi = rnorm(d, 0, sqrt(M))
  th_old = th
  
  log_p_old = log_posterior(X, y, th) - .5*sum(M_inv * phi^2)
  
  phi = phi + .5 * epsilon * gradient_theta(X, y, th)
  
  for (l in 1:L) {
    th  = th + epsilon*M_inv*phi
    phi = phi + ifelse(l == L, .5, 1) * epsilon * gradient_theta(X, y, th)
  }
  
  # here we get into standard MCMC stuff, jump or not based on a draw from a
  # proposal distribution
  phi = -phi
  log_p_star = log_posterior(X, y, th) - .5*sum(M_inv * phi^2)    
  r = exp(log_p_star - log_p_old)
  
  if (is.nan(r)) r = 0
  
  p_jump = min(r, 1)
  
  if (runif(1) < p_jump) {
    th_new = th
  }
  else {
    th_new = th_old
  }
  
  # returns estimates and acceptance rate
  list(th = th_new, p_jump = p_jump)  
}

##---- 5th part
hmc_run <- function(starts, iter, warmup, epsilon_0, L_0, M, X, y) {
  # # Args: 
  # starts:  starting values
  # iter: total number of simulations for each chain (note chain is based on the dimension of starts)
  # warmup: determines which of the initial iterations will be ignored for inference purposes
  # epsilon0: the baseline stepsize
  # L0: the baseline number of leapfrog steps
  # M: is the mass vector
  chains = nrow(starts)
  d = ncol(starts)
  sims = array(NA, 
               c(iter, chains, d), 
               dimnames = list(NULL, NULL, colnames(starts)))
  p_jump = matrix(NA, iter, chains)
  
  for (j in 1:chains) {
    th = starts[j,]
    
    for (t in 1:iter) {
      if (t%%50==0){
        print(t)
      }
      epsilon = runif(1, 0, 2*epsilon_0)
      L    = ceiling(2*L_0*runif(1))
      
      temp = hmc_iteration(X, y, th, epsilon, L, M)
      
      p_jump[t,j] = temp$p_jump
      sims[t,j,]  = temp$th
      
      th = temp$th
    }
  }
  
  # acceptance rate
  acc = round(colMeans(p_jump[(warmup + 1):iter,]), 3)  
  
  message('Avg acceptance probability for each chain: ', 
          paste0(acc[1],', ',acc[2]), '\n') 
  
  list(sims = sims, p_jump = p_jump)
}

##---- 6th
# Starting values and mcmc settings
parnames = c(paste0('beta[', 1:4, ']'), 'sigma')
d = length(parnames)

chains = 2

theta_start = t(replicate(chains, c(runif(d-1, -1, 1), 1)))
colnames(theta_start) = parnames

nsim = 1000
wu   = 500

stepsize = .08
nLeap = 10
vars  = rep(1, 5)
mass_vector = 1 / vars


##---- 7th part
# Run the model
fit_hmc = hmc_run(
  starts    = theta_start,
  iter      = nsim,
  warmup    = wu,
  epsilon_0 = stepsize,
  L_0 = nLeap,
  M   = mass_vector,
  X   = X,
  y   = y
)
# str(fit_hmc, 1)

##---- 8th

library(coda)

theta = as.mcmc.list(list(as.mcmc(fit_hmc$sims[(wu+1):nsim, 1,]), 
                          as.mcmc(fit_hmc$sims[(wu+1):nsim, 2,])))

# summary(theta)
fit_summary =  summary(theta)$statistics[,'Mean']

beta_est  = fit_summary[1:4]
sigma_est = fit_summary[5]

log_posterior(X, y, fit_summary)

# Instead we can use rstan’s monitor function on fit_hmc$sims to produce typical Stan output.
library(rstan)
monitor(fit_hmc$sims)