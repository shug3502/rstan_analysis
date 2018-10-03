## trying mixture model with synthetic data to check for identifiability
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

identifier = 'syntheticv001'
run_mcmc=TRUE
nSamples = 10
nTest=1
nTestOE=0

source('extract_times_and_scaling.R')
times = extract_times_and_scaling(nSamples,nTest,nTestOE)
#############################################################

expose_stan_functions('model_comparison_at_stst3.stan')
#suppose we observe data from one of the mixture components and feed this as data to the mixture model
#can we identify the correct mixture component?
nu_real = 0.95
phi_real = 0.35
xi_real = 0.1
rc_index_real = 2

raw_data = get_k2(nu_real,get_RC_from_dict(rc_index_real))
data = (raw_data * c(1,rep(1/phi_real,15))) %>% rep(.,nSamples) %>% matrix(.,ncol=16,byrow = TRUE) + matrix(abs(xi_real*rnorm(16*nSamples)),ncol=16)
plot(data[1,])
##############################################################


if (run_mcmc) {
  stan_file <- 'model_comparison_at_stst3.stan'
  estimates <- stan(file = stan_file,
                    data = list (
                      T1 = nSamples,
                      T2 = nTest+nTestOE,
                      y_obs = data
                    ),
                    seed = 42,
                    chains = 4,
                    warmup = 1000,
                    iter = 2000,
                    init = function() list(nu=0.9, xi=0.4, phi=0.3)                    
  )
  
  tryCatch({
    saveRDS(estimates, file = paste('fits/model_comparison',identifier,'.rds',sep='')) # to load use readRDS("fit.rds")
  }, warning = function(war) {  
    # warning handler picks up where error was generated
    print(paste("WARNING:  ",war))
  }, error = function(err) {  
    # error handler picks up where error was generated
    print(paste("NOT SAVED ERROR:  ",err))
  }) # END tryCatch
} else {
  estimates = readRDS(paste('fits/model_comparison',identifier,'.rds',sep=''))
}
print(estimates, pars = c('nu','phi','xi'))

