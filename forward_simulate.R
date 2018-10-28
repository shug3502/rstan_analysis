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
source('get_blocked_indices.R')
expose_stan_functions('mrna_transport_with_blocking.stan')
###############
identifier_simple = 'v470_with_outliers_M0_simple' #'v431WT_simple' #load fitted posterior parameters from this model
identifier_DD = 'v470_with_outliers_M8_density_dependent_4a'
nTestOE = 14
num_draws = 1000
models_to_sim = paste('M',seq(from=0,to=14)[-6],sep='')
models_with_blocking = c('M1','M6')
models_with_density_dependence = c('M7','M8','M9')
times = extract_times_and_scaling(1,1,nTestOE)
multipliers <- list(M0=c(2,2),M1=c(2,2),M2=c(2,1),M3=c(4,2),M4=c(4,1),M6=c(4,4),M7=c(2,2),M8=c(4,4),M9=c(4,4),M13=c(4,4),M14=c(2.5,2.5))
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

#informal test of log likelihood
th = c(0.2, 10.0,  0.9, 1.0,  0.3,  1.3)
forward_simulate_wrapper(0.2,10,0.9,0.01,0.3,1.3,'M0')
get_log_lik(overexpression_data,nTestOE,th,producers_list[['M0']],blocked_list[['M0']],times)

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

post_pred_simple <- draws_simple %>%
  mutate(modelID = list(models_to_sim[!(models_to_sim %in% models_with_density_dependence)])) %>%
  unnest() %>%
  mutate(time=list(times$ts4),
         y_pred_oe = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
                                 forward_simulate_wrapper)) %>%
  unnest(y_pred_oe,time,.preserve=modelID) %>%
  mutate(cellID = list(seq_len(16))) %>%
  unnest(.preserve=modelID) %>%
  unnest() %>%
  group_by(time,cellID,modelID) %>%
  mutate(y_lower=quantile(y_pred_oe,0.025),
         y_upper=quantile(y_pred_oe,0.975),
         y_median=quantile(y_pred_oe,0.5)) %>%
  ungroup()

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

###################
expose_stan_functions('mrna_transport_density_dependent_with_blocking.stan') #should replace previous functions and work with previous wrappers
draws_DD <- estimates_DD %>%
  spread_draws(a,b,nu,beta,phi,sigma) %>% 
  filter(.iteration<=num_draws) %>%
  select(b,a,nu,beta,phi,sigma)

post_pred_DD <- draws_DD %>%
  mutate(modelID = list(models_to_sim[models_to_sim %in% models_with_density_dependence])) %>%
  unnest() %>%
  mutate(time=list(times$ts4),
         y_pred_oe = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
                                 forward_simulate_wrapper)) %>%
  unnest(y_pred_oe,time,.preserve=modelID) %>%
  mutate(cellID = list(seq_len(16))) %>%
  unnest(.preserve=modelID) %>%
  unnest() %>%
  group_by(time,cellID,modelID) %>%
  mutate(y_lower=quantile(y_pred_oe,0.025),
         y_upper=quantile(y_pred_oe,0.975),
         y_median=quantile(y_pred_oe,0.5)) %>%
  ungroup()

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

full_log_lik_df <- full_join(log_lik_simple, log_lik_DD)

log_lik_df2 <- full_log_lik_df %>%
  filter(cellID!=1) %>%
  group_by(modelID,b,a,nu,phi,sigma) %>%
  summarise(ll_sum = sum(ll))

######################

#add observed overexpression data
xdata <- data.frame(rna = as.vector(overexpression_data),cellID = as.vector(matrix(rep(1:16,nTestOE),nrow=nTestOE,byrow=TRUE)),time = rep(times$ts4,16))
observed_and_predictions <- full_join(full_join(post_pred_simple, post_pred_DD), xdata)

observed_and_predictions %>%
  filter(round(time)==28) %>%
  ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
  geom_ribbon(alpha=0.4) +
  geom_line() +
  geom_point(aes(x=cellID,y=rna)) +
  facet_wrap(~factor(modelID))

observed_and_predictions %>%
  filter(modelID=='M6') %>%
  ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
  geom_ribbon(alpha=0.4) +
  geom_line() +
  geom_point(aes(x=cellID,y=rna)) +
  facet_wrap(~factor(time),scales='free_y')  

g <- log_lik_df2 %>%
  ggplot(aes(ll_sum,color=modelID)) +
  geom_density() + 
  scale_x_continuous(limits=c(-1500,-800))
  print(g)
  ggsave(paste('plots/compare_production_models_fwd_sim_',identifier_simple,'.eps',sep=''),device=cairo_ps)
  
  log_lik_df2 %>%
    filter(as.numeric(substr(modelID,2,3))>10 | modelID=='M6' | modelID=='M3') %>%
    ggplot(aes(ll_sum,color=modelID)) +
    geom_density() + 
    scale_x_continuous(limits=c(-1100,-850))
  
  ###########################
  
  analyse_given_model <- function(full_log_lik,nColumns,mID='M0'){
    #mID is a string containing the model ID
    library(loo)
    log_lik_matrix = matrix(full_log_lik %>% filter(modelID==mID) %>% .$ll_sum,ncol=nColumns,byrow=FALSE)
    loo1 <- loo(log_lik_matrix)
    return(loo1) 
  }

  #look at LOO by cellID
  log_lik_by_cellID <- full_log_lik_df %>%
    # filter(cellID!=1) %>% 
    group_by(cellID,modelID,a,b,phi,sigma,nu) %>%
    summarise(ll_sum=sum(ll))
  aM_cellID_list <- purrr::map(models_to_sim, function(x) analyse_given_model(log_lik_by_cellID,16,x))
  purrr::map(aM_cellID_list,plot)

  #look at LOO by egg chamber or time pt
  log_lik_by_time <- full_log_lik_df %>%
    # filter(cellID!=1) %>%     
    group_by(time,modelID,a,b,phi,sigma,nu) %>%
    summarise(ll_sum=sum(ll))
  aM_time_list <- purrr::map(models_to_sim, function(x) analyse_given_model(log_lik_by_time,nTestOE,x))
  purrr::map(aM_time_list,plot)
  