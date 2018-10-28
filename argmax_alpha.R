#maxmise the expected loglikelihood value of a model
#given a multiplier alpha time wild type production
#################

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
source('get_blocked_indices.R')
expose_stan_functions('mrna_transport_with_blocking.stan')
##############

calculate_expected_ll <- function(alpha,
                                  TestOE = 9,
                                  num_draws = 100,
                                  identifier_simple = 'v470_with_outliers_M0_simple',
                                  identifier_DD = 'v470_with_outliers_M8_density_dependent_4a',
                                  models_to_sim = paste('Malpha'),
                                  with_blocking = FALSE,
                                  with_density_dependence = FALSE
                                  ){
if (with_blocking){
  models_with_blocking = c('M1','M6','Malpha')
} else {
  models_with_blocking = c('M1','M6')
}
if (with_density_dependence){
  models_with_density_dependence = c('M7','M8','M9','Malpha')
} else {
  models_with_density_dependence = c('M7','M8','M9')
}
times = extract_times_and_scaling(1,1,nTestOE)
multipliers <- list(Malpha=c(alpha,alpha),M0=c(2,2),M1=c(2,2),M2=c(2,1),M3=c(4,2),M4=c(4,1),M6=c(4,4),M7=c(2,2),M8=c(4,4),M9=c(4,4),M13=c(4,4),M14=c(2.5,2.5))
blocked_matrix = get_blocked_indices(nTestOE)[times$sort_indices4,] #matrix of which RCs are blocked
overexpression_data = matrix(as.numeric(read.csv('data/exp_data_overexpression.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
overexpression_data = overexpression_data[times$sort_indices4,] #need to sort time series and correspondingly reorder rows

producers_list <- list()
blocked_list <- list()
for (modelID in models_to_sim) {
  if (modelID %in% names(multipliers)){
    producers_list[modelID] = list(get_producers(nTestOE,multiplier=multipliers[[modelID]]))
  } else if (modelID=='M10') {
    producers_list[modelID] = list(get_adjusted_producers(nTestOE))
  } else if (modelID=='M11') {
    producers_list[modelID] = list(estimate_adjusted_producers(nTestOE))
  } else if (modelID=='M12') {
    producers_list[modelID] = list(2*estimate_adjusted_producers(nTestOE))
  } else {
    stop(paste('Oops, do not have a multiplier for model: ',modelID, sep=''))
  }
  if (modelID %in% models_with_blocking){
    blocked_list[modelID] = list(blocked_matrix)
  } else {
    blocked_list[modelID] = list(matrix(rep(1,3*nTestOE),ncol=3))
  }
}

#######################

forward_simulate <- function(nTestOE,theta,producers,blocked,times){
  #nTestOE: int, number of time points of OE data
  #theta: numeric vector, parameter values to evaluate at, for example drawn from posterior distribution (b,a,nu,phi,sigma)
  #producers: numeric matrix, size nTestOE x 16, each row contains info on production in each nurse cell  
  #times: list, produced from extract_times_and_scaling, contains time series to evaluate OE data at
  ###########
  blocked_int_list = lapply(1:nTestOE, FUN = function(i) blocked[i,]) 
  y_pred_oe <- forward_simulate_OE_rng(nTestOE,times$ts4,rep(0,16),theta,producers,blocked_int_list)
  return(y_pred_oe) 
}

forward_simulate_wrapper <- function(b,a,nu,beta,phi,sigma,modelID='M0'){
  #model ID indicates model used for producers
  #M0 is (2a,2a) production
  #M2 is (2a,a)
  #M3 is (4a,2a)
  #M4 is (4a,a)
  #M10 is data driven without missing data model
  #M11 is estimated from data with missing data model
  ###############
  th = c(b,a,nu,beta,phi,sigma)
  producers = producers_list[[modelID]]
  blocked = blocked_list[[modelID]]
  s = forward_simulate(nTestOE,th,producers,blocked,times)
  return(s)
}

get_log_lik <- function(overexpression_data,nTestOE,theta,producers,blocked,times){
  overexpression_int_list = lapply(1:nTestOE, FUN = function(i) overexpression_data[i,]) #Rcpp is specific about the input format of arguments, so haveto make things into lists
  blocked_int_list = lapply(1:nTestOE, FUN = function(i) blocked[i,])   #similarly
  ll <- my_log_lik(overexpression_int_list,nTestOE,times$ts4,rep(0,16),theta,producers,blocked_int_list)
  return(ll) 
}

get_log_lik_wrapper <- function(b,a,nu,beta,phi,sigma,modelID='M0'){
  th = c(b,a,nu,beta,phi,sigma)
  producers = producers_list[[modelID]]
  blocked = blocked_list[[modelID]]
  ll = get_log_lik(overexpression_data,nTestOE,th,producers,blocked,times)
  #easier for getting into right form in dataframe later (via unnest) to put back into list
  return(lapply(1:nTestOE, FUN = function(i) ll[i,]))
}

check_and_load_path <- function(fit_path){
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
  return(estimates)
}

#get results from fitting simple model
cat(paste('\nloading fitted simple model parameters from model id: ',identifier_simple,'\n',sep=''))
path_simple = paste('fits/mrna_transport_estimates',identifier_simple,'.rds',sep='')
estimates_simple <- check_and_load_path(path_simple)
cat(paste('\nloading fitted density dependent model parameters from model id: ',identifier_DD,'\n',sep=''))
path_DD = paste('fits/mrna_transport_estimates',identifier_DD,'.rds',sep='')
estimates_DD <- check_and_load_path(path_DD)
#################

draws_simple <- estimates_simple %>%
  spread_draws(a,b,nu,phi,sigma) %>% 
  filter(.iteration<=num_draws) %>%
  select(b,a,nu,phi,sigma) %>%
  mutate(beta = -1) 

log_lik_simple <- draws_simple %>%
  mutate(modelID = list(models_to_sim[!(models_to_sim %in% models_with_density_dependence)])) %>%
  unnest() %>%
  mutate(time=list(times$ts4),
         ll = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
                          get_log_lik_wrapper)) %>%
  unnest(ll,time,.preserve=modelID) %>%
  mutate(cellID = list(seq_len(16))) %>%
  unnest(ll,cellID,.preserve=modelID) %>%
  unnest()

##################

expose_stan_functions('mrna_transport_density_dependent_with_blocking.stan') #should replace previous functions and work with previous wrappers
draws_DD <- estimates_DD %>%
  spread_draws(a,b,nu,beta,phi,sigma) %>% 
  filter(.iteration<=num_draws) %>%
  select(b,a,nu,beta,phi,sigma)

log_lik_DD <- draws_DD %>%
  mutate(modelID = list(models_to_sim[models_to_sim %in% models_with_density_dependence])) %>%
  unnest() %>%
  mutate(time=list(times$ts4),
         ll = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
                          get_log_lik_wrapper)) %>%
  unnest(ll,time,.preserve=modelID) %>%
  mutate(cellID = list(seq_len(16))) %>%
  unnest(ll,cellID,.preserve=modelID) %>%
  unnest()

#####################

expected_log_lik <- full_join(log_lik_simple, log_lik_DD) %>%
  group_by(modelID,b,a,nu,phi,sigma,beta) %>%
  summarise(sum_ll=sum(ll)) %>%
  group_by(modelID) %>%
  summarise(ell = mean(sum_ll)) %>%
  arrange(desc(ell))

}
