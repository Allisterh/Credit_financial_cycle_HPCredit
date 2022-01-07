clear, clc
working_dir = '/Users/namnguyen/Documents/GitHub/HPCredit/Regression/Bayesian_UC_AR2/HMC/'
cd(working_dir)
%************ Data generation
N = 250
K = 3

covariates = normrnd(0,1, [N,K]);

%** Create the model matrix with intercept
X=[ones(N,1) covariates];


% create a normally distributed variable that is a function of the covariates
coefs = [5, .2, -1.5, .9]';
sigma = 2;
mu = X * coefs;
y  = normrnd(mu, sigma, [N,1]);


% Starting values and mcmc settings
% parnames = c(paste0('beta[', 1:4, ']'), 'sigma')
d = 5;

chains = 2;

theta_start = [unifrnd(-1,1, [2,d-1]), repmat(1,[2,1])];

nsim = 1000;
wu   = 500;

stepsize = .08;
nLeap = 10;
vars  = repmat(1, 1, 5);
mass_vector = 1./vars; 


%%---- 7th part
% Run the model
fit_hmc = hmc_run(  theta_start, nsim, wu,  stepsize,  nLeap,  mass_vector,  X,  y);

phi= mvnrnd(zeros(5,1), sqrt(mass_vector))

% str(fit_hmc, 1)

%---- 8th

library(coda)

theta = as.mcmc.list(list(as.mcmc(fit_hmc$sims[(wu+1):nsim, 1,]), 
                          as.mcmc(fit_hmc$sims[(wu+1):nsim, 2,])))

% summary(theta)
fit_summary =  summary(theta)$statistics[,'Mean']

beta_est  = fit_summary[1:4]
sigma_est = fit_summary[5]

log_posterior(X, y, fit_summary)

% Instead we can use rstan’s monitor function on fit_hmc$sims to produce typical Stan output.
library(rstan)
monitor(fit_hmc$sims)