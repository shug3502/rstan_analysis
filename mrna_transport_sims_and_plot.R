#setwd('~/Documents/FISH_data/rstan_analysis')
mrna_transport_sims_and_plot <- function(identifier='v003',nSamples=15,nTest=5,th=c(2,2),sig=10^-9,phi=0.289,nu=0.95){
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# identifier = 'v003' #run identifier
# nSamples = 15 #how many egg chambers segmented
# nTest = 5

#############################################################
#Fit linear model to log(length) of egg chambers to get time of development.
#Could use area or volume of egg chamber. Compare to Jia 2016 Sci Reports.
#Area may be better as in some examples egg chambers are quite squashed.
#To get lengths read length.txt file in each folder
source('extract_times_and_scaling.R')
times = extract_times_and_scaling(nSamples,nTest)

#############################################################
m0 = c(0, rep(1,15)) #initial condition
#th = c(2,2)
#sig = 10^-9
#phi = 0.289

t0 = times$t0 #0.0
ts = times$ts1 
#nu = 0.95

#############################################################
mc <- stan_model('model_comparison5.stan')
expose_stan_functions(mc)
B = construct_matrix(nu) %>% as.vector
#############################################################

  #sample from the model to get fake data
  #print('using fake generated data')
  samples <- stan(file = 'mrna_transport6.stan',
                  data = list (
                    T  = nSamples,
                    y0 = m0,
                    t0 = t0,
                    ts = ts,
                    theta = array(th, dim = 2),
                    sigma = sig,
                    phi = phi,
                    B = B
                  ),
                  algorithm="Fixed_param",
                  seed = 42,
                  chains = 1,
                  iter =100, 
                  refresh = -1
  )
  
source('post_pred_plot.R')
p1 <- post_pred_plot(raw_data=NA,ts,nSamples,'y_hat',samples,identifier,title_stem='plots/sims')
print(p1)
return(p1)
}