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
models_to_sim = list('M0','M2','M3','M4','M10','M11','M12','M13','M14')
num_draws = 100
multipliers <- list(M0=c(2,2),M2=c(2,1),M3=c(4,2),M4=c(4,1),M13=c(4,4),M14=c(2.5,2.5))
identifier = 'v431WT_simple' #'v470_with_outliers_M0_simple' #load fitted posterior parameters from this model
times = extract_times_and_scaling(1,1,nTestOE)
overexpression_data = matrix(as.numeric(read.csv('data/exp_data_overexpression.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
overexpression_data = overexpression_data[times$sort_indices4,] #need to sort time series and correspondingly reorder rows

producers_list <- list()
for (modelID in models_to_sim) {
  if (modelID %in% names(multipliers)){
    producers_list[modelID] = list(get_producers(nTestOE,multiplier=multipliers[[modelID]]))
  } else if (modelID=='M10') {
    producers_list[modelID] = list(get_adjusted_producers(nTestOE))
  } else if (modelID=='M11') {
    producers_list[modelID] = list(estimate_adjusted_producers(nTestOE))
  } else if (modelID=='M12') {
    producers_list[modelID] = list(2*estimate_adjusted_producers(nTestOE))
}
}

#######################

forward_simulate <- function(nTestOE,theta,producers,times){
  #nTestOE: int, number of time points of OE data
  #theta: numeric vector, parameter values to evaluate at, for example drawn from posterior distribution (b,a,nu,phi,sigma)
  #producers: numeric matrix, size nTestOE x 16, each row contains info on production in each nurse cell  
  #times: list, produced from extract_times_and_scaling, contains time series to evaluate OE data at
  ###########
  y_pred_oe <- forward_simulate_OE_rng(nTestOE,times$ts4,rep(0,16),theta,producers)
  return(y_pred_oe) 
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
  s = forward_simulate(nTestOE,th,producers,times)
  return(s)
  }

get_log_lik <- function(overexpression_data,nTestOE,theta,producers,times){
  overexpression_data_list = lapply(1:nTestOE, FUN = function(i) overexpression_data[i,])
  ll <- my_log_lik(overexpression_data_list,nTestOE,times$ts4,rep(0,16),theta,producers)
  return(ll) 
}

get_log_lik_wrapper <- function(b,a,nu,phi,sigma,modelID='M0'){
  th = c(b,a,nu,phi,sigma)
  producers = producers_list[[modelID]]
  s = get_log_lik(overexpression_data,nTestOE,th,producers,times)
  return(sum(s))
}

#test log likelihood
th = c(0.2, 10.0,  0.9,  0.3,  1.3)
get_log_lik(overexpression_data,nTestOE,th,producers_list[['M0']],times)


#get results from fitting model
cat(paste('\nloading fitted model parameters from model id: ',identifier,'\n',sep=''))
fit_path = paste('fits/mrna_transport_estimates',identifier,'.rds',sep='')
if (file.exists(fit_path)){
  estimates <- readRDS(fit_path)
} else {
  #try on scratch
  fit_path = paste('/scratch/harrison/FISH_data/old_fits/mrna_transport_estimates',identifier,'.rds',sep='')
  if (file.exists(fit_path)){
    estimates <- readRDS(fit_path)
  } else {
    stop('file does not exist, please run analysis or check a different computer')
  }
}

#################

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

log_lik_df <- y_pred0pt0 %>%
  mutate(modelID = list(models_to_sim)) %>%
  unnest(.preserve=theta) %>%
  mutate(time=list(times$ts4),
         ll = purrr::pmap(list(b,a,nu,phi,sigma,modelID),
                                 get_log_lik_wrapper)) %>%
  unnest(time,.preserve=c(ll,modelID)) %>%
  unnest(ll,.preserve = modelID) %>%
  unnest() %>%
  mutate(cellID = list(seq_len(16))) %>%
  unnest()

log_lik_df2 <- y_pred0pt0 %>%
  mutate(modelID = list(models_to_sim)) %>%
  unnest(.preserve=theta) %>%
  mutate(ll = purrr::pmap(list(b,a,nu,phi,sigma,modelID),
                          get_log_lik_wrapper)) %>%
  unnest(modelID,ll)

post_pred_df %<>%
  unnest() %>%
  group_by(time,cellID,modelID) %>%
  mutate(y_lower=quantile(y_pred_oe,0.025),
         y_upper=quantile(y_pred_oe,0.975),
         y_median=quantile(y_pred_oe,0.5)) %>%
  ungroup()

#add observed overexpression data
xdata <- data.frame(rna = as.vector(overexpression_data),cellID = as.vector(matrix(rep(1:16,nTestOE),nrow=nTestOE,byrow=TRUE)),time = rep(times$ts4,16))

observed_and_predictions <- full_join(post_pred_df, xdata)

observed_and_predictions %>%
  filter(round(time)==28) %>%
  ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
  geom_ribbon(alpha=0.4) +
  geom_line() +
  geom_point(aes(x=cellID,y=rna)) +
  facet_wrap(~factor(modelID))

observed_and_predictions %>%
  filter(modelID=='M12') %>%
  ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
  geom_ribbon(alpha=0.4) +
  geom_line() +
  geom_point(aes(x=cellID,y=rna)) +
  facet_wrap(~factor(time))  

g <- log_lik_df2 %>%
  unnest(ll) %>%
  ggplot(aes(ll,color=modelID)) +
  geom_density() + 
  scale_x_continuous(limits=c(-1800,-800))
  print(g)
  ggsave(paste('plots/compare_production_models_fwd_sim_',identifier,'.eps'),device=cairo_ps)
  
  log_lik_df2 %>%
    filter(as.numeric(substr(modelID,2,3))>10 | modelID=='M3') %>%
    unnest(ll) %>%
    ggplot(aes(ll,color=modelID)) +
    geom_density()
  