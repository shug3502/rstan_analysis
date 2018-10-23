## Inhomogeneous producers for overexpressor seems a plausible model, but implementing this in a sensible way is hard
## In particular choice of distribution of RNA producers across cells and total RNA production int eh egg chamber in relation to wild type
##
## possible choices are uniform distribution of producers across cells, or data driven versions
## 2*15*a total production whre wild type production is 15*a
##
###############
## Here we implement a function(s) that takes an array of model parameters and simulates from the forward model to give posterior predictive samples

###########
#setup
library(rstan)
library(dplyr)
library(tidyr)
library(magrittr)
library(tidybayes)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('get_producers.R')
source('get_adjusted_producers.R')
source('estimate_adjusted_producers.R')
source('extract_times_and_scaling.R')
expose_stan_functions('M0.stan')
###############
nTestOE = 9
models_to_sim = list('M0','M2','M3','M4','M10','M11')
num_draws = 50
multipliers <- list(M0=c(2,2),M2=c(2,1),M3=c(4,2),M4=c(4,1))
identifier = 'v431WT_simple' #load fitted posterior parameters from this model
times = extract_times_and_scaling(1,1,nTestOE)

producers_list <- list()
for (modelID in models_to_sim) {
  if (modelID %in% names(multipliers)){
    producers_list[modelID] = list(get_producers(nTestOE,multiplier=multipliers[[modelID]]))
  } else if (modelID=='M10') {
    producers_list[modelID] = list(get_adjusted_producers(nTestOE))
  } else if (modelID=='M11') {
    producers_list[modelID] = list(estimate_adjusted_producers(nTestOE))
  }
}

#######################

forward_simulate_M0 <- function(nTestOE,theta,producers,times){
  #nTestOE: int, number of time points of OE data
  #theta: numeric vector, parameter values to evaluate at, for example drawn from posterior distribution (b,a,nu,phi,sigma)
  #producers: numeric matrix, size nTestOE x 16, each row contains info on production in each nurse cell  
  #times: list, produced from extract_times_and_scaling, contains time series to evaluate OE data at
  ###########
  y_pred_oe <- forward_simulate_OE_rng(nTestOE,times$ts4,rep(0,16),theta,producers)
  return(y_pred_oe)  #(matrix(unlist(y_pred_oe),ncol=16,byrow=TRUE))
}

forward_simulate_wrapper <- function(b,a,nu,phi,sigma,modelID='M0'){
  #model ID indicates model used for producers
  #M0 is (2a,2a) production
  #M2 is (2a,a)
  #M3 is (4a,2a)
  #M4 is (4a,a)
  #M10 is data driven without missing data model
  #M11 is estimated from data with missing data model
  ###############
  th = c(b,a,nu,phi,sigma)
  producers = producers_list[[modelID]]
  s = forward_simulate_M0(nTestOE,th,producers,times)
  return(s)
  }

#get results from fitting model
cat(paste('\nloading fitted model parameters from model id: ',identifier,'\n',sep=''))
estimates <- readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep=''))
# e = rstan::extract(estimates,pars=c('theta','phi','sigma'))
# theta = cbind(e[['theta']],e[['phi']],e[['sigma']]) #probably easier way to do this

##################
# #Now try different models
# #M0
# prod0 <- matrix(rep(2,16*nTestOE),ncol=16)
# prod0[,1] = 0
# #another way to do this by having all producers in a list and iterating over
# y_pred0 <- apply(theta[seq_len(100),],1,function(x) forward_simulate_M0(nTestOE,x,prod0,times))
# 
# #M2
# source('get_producers.R')
# prod2 <- get_producers(nTestOE,multiplier=c(2,1))
# y_pred2 <- apply(theta[seq_len(100),],1,function(x) forward_simulate_M0(nTestOE,x,prod2,times))
# 
# #M3
# prod3 <- get_producers(nTestOE,multiplier=c(4,2))
# y_pred3 <- apply(theta[seq_len(100),],1,function(x) forward_simulate_M0(nTestOE,x,prod3,times))
# 
# 
# #M4
# prod4 <- get_producers(nTestOE,multiplier=c(4,1))
# y_pred4 <- apply(theta[seq_len(100),],1,function(x) forward_simulate_M0(nTestOE,x,prod4,times))
# 
# 
# #M10
# source('get_adjusted_producers.R')
# prod10 <- get_adjusted_producers(nTestOE)
# y_pred10 <- apply(theta[seq_len(100),],1,function(x) forward_simulate_M0(nTestOE,x,prod10,times))
# 
# 
# #M11
# source('estimate_adjusted_producers.R')
# prod11 = estimate_adjusted_producers(nTestOE)
# y_pred11 <- apply(theta[seq_len(100),],1,function(x) forward_simulate_M0(nTestOE,x,prod11,times))
# 
# 
# #want a data frame with the following fields
# #cellID, time, modelID, theta, predicted_RNA, production
# 
# 
# 
# 
# #what do we want to do with the output now?
# #we have models and 



############
#leave tidy version for now

#same kind of thing using tidybayes

y_pred0pt0 <- estimates %>%
  spread_draws(a,b,nu,phi,sigma) %>% 
  filter(.iteration<=num_draws) %>%
  select(b,a,nu,phi,sigma) %>%
  mutate(theta = purrr::pmap(list(b,a,nu,phi,sigma),function(b,a,nu,phi,sigma) c(b,a,nu,phi,sigma)))
post_pred_df <- y_pred0pt0 %>%
  mutate(modelID = list(models_to_sim)) %>%
  unnest(.preserve=theta) %>%
  mutate(time=list(times$ts4),
         y_pred_oe = purrr::pmap(list(b,a,nu,phi,sigma,modelID),
                                 forward_simulate_wrapper)) %>%
  unnest(y_pred_oe,time,.preserve=modelID) %>%
  mutate(cellID = list(seq_len(16))) %>%
  unnest(.preserve=modelID)

post_pred_df %<>%
  unnest() %>%
  group_by(time,cellID,modelID) %>%
  mutate(y_lower=quantile(y_pred_oe,0.025),
         y_upper=quantile(y_pred_oe,0.975),
         y_median=quantile(y_pred_oe,0.5)) %>%
  ungroup()

post_pred_df %>%
  filter(round(time)==19) %>%
  ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
  geom_ribbon(alpha=0.4) +
  geom_line() +
  facet_wrap(~factor(modelID))

post_pred_df %>%
  filter(modelID=='M0') %>%
  ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
  geom_ribbon(alpha=0.4) +
  geom_line() +
  facet_wrap(~factor(time))  

  