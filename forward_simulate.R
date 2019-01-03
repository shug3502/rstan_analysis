## Inhomogeneous producers for overexpressor seems a plausible model, but implementing this in a sensible way is hard
## In particular choice of distribution of RNA producers across cells and total RNA production int eh egg chamber in relation to wild type
##
## possible choices are uniform distribution of producers across cells, or data driven versions
## 2*15*a total production whre wild type production is 15*a
##
###############
## Here we implement a function(s) that takes an array of model parameters and simulates from the forward model to give posterior predictive samples
##This is set up to load and deal with the model with inhomogeneous overdispersion parameter sigma
## inhomogeneous sigma
###########
#setup
library(rstan)
library(dplyr)
library(tidyr)
library(magrittr)
library(tidybayes)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('get_producers.R')
source('get_adjusted_producers.R')
source('estimate_adjusted_producers.R')
source('extract_times_and_scaling.R')
source('get_blocked_indices.R')
expose_stan_functions('mrna_transport_with_blocking.stan')
###############
identifier_simple = 'v510_minimalM0_simple' #'v470_with_outliers_M0_simple' #'v431WT_simple' #load fitted posterior parameters from this model
identifier_DD = 'v510_minimalM4_density_dependent'  #'v470_with_outliers_M8_density_dependent_4a'
nTestOE = 9
num_draws = 100
gamma = 2
models_to_sim = paste('M',seq(from=0,to=7),sep='')
models_with_blocking = c('M1','M5','M6','M7')
models_with_density_dependence = c('M2','M4','M6','M7')
models_with_data_driven_production = c('M3','M4','M5','M7')
times = extract_times_and_scaling(1,1,nTestOE)
multipliers <- list(M0=gamma,M1=gamma,M2=gamma,M6=gamma)
blocked_matrix = get_blocked_indices(nTestOE)[times$sort_indices4,] #matrix of which RCs are blocked
overexpression_data = matrix(as.numeric(read.csv('data/exp_data_overexpression.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
overexpression_data = overexpression_data[times$sort_indices4,] #need to sort time series and correspondingly reorder rows

producers_list <- list()
blocked_list <- list()
for (modelID in models_to_sim) {
  if (modelID %in% names(multipliers)){
    producers_list[modelID] = list(get_producers(nTestOE,multiplier=multipliers[[modelID]]))
  } else if (modelID %in% models_with_data_driven_production) {
    producers_list[modelID] = list(estimate_adjusted_producers(nTestOE,gamma=gamma))
  } else if (modelID=='M10') {
    producers_list[modelID] = list(get_adjusted_producers(nTestOE))
  } else if (modelID=='M12') {
    producers_list[modelID] = list(estimate_adjusted_producers(nTestOE,gamma=4))
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

forward_simulate <- function(nTestOE,theta,sigma,producers,blocked,times){
  #nTestOE: int, number of time points of OE data
  #theta: numeric vector, parameter values to evaluate at, for example drawn from posterior distribution (b,a,nu,phi,sigma)
  #producers: numeric matrix, size nTestOE x 16, each row contains info on production in each nurse cell  
  #times: list, produced from extract_times_and_scaling, contains time series to evaluate OE data at
  ###########
  blocked_int_list = lapply(1:nTestOE, FUN = function(i) blocked[i,]) 
  y_pred_oe <- forward_simulate_OE_rng(nTestOE,times$ts4,rep(0,16),theta,sigma,producers,blocked_int_list)
  return(y_pred_oe) 
}

forward_simulate_wrapper <- function(b,a,nu,beta,phi,
                                     sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,
                                     sigma9,sigma10,sigma11,sigma12,sigma13,sigma14,sigma15,sigma16,
                                     modelID='M0'){
  #model ID indicates model used for producers
  ###############
  th = c(b,a,nu,beta,phi)
  sigma = c(sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,
            sigma9,sigma10,sigma11,sigma12,sigma13,sigma14,sigma15,sigma16)
  producers = producers_list[[modelID]]
  blocked = blocked_list[[modelID]]
  s = forward_simulate(nTestOE,th,sigma,producers,blocked,times)
  return(s)
  }

get_log_lik <- function(overexpression_data,nTestOE,theta,sigma,producers,blocked,times){
  overexpression_int_list = lapply(1:nTestOE, FUN = function(i) overexpression_data[i,]) #Rcpp is specific about the input format of arguments, so haveto make things into lists
  blocked_int_list = lapply(1:nTestOE, FUN = function(i) blocked[i,])   #similarly
  ll <- my_log_lik(overexpression_int_list,nTestOE,times$ts4,rep(0,16),theta,sigma,producers,blocked_int_list)
  return(ll) 
}

alternative_log_lik_wrapper <- function(b,a,nu,beta,phi,
                                        sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,
                                        sigma9,sigma10,sigma11,sigma12,sigma13,sigma14,sigma15,sigma16,
                                        modelID='M0',drop_egg_ind=FALSE,drop_cell_ind=FALSE){
  th = c(b,a,nu,beta,phi)
  # I'm guessing theres a better way than taking each as a separate column from the dataframe
  sigma = c(sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,
            sigma9,sigma10,sigma11,sigma12,sigma13,sigma14,sigma15,sigma16)
  producers = producers_list[[modelID]]
  blocked = blocked_list[[modelID]]
  ll = get_log_lik(overexpression_data,nTestOE,th,sigma,producers,blocked,times)
  if (drop_egg_ind){ #allows us to look at the effect of leaving out an egg chamber or cell 
    ll = ll[-drop_egg_ind,]
  }
  if (drop_cell_ind){
    ll = ll[,-drop_cell_ind]
  }
  #easier for getting into right form in dataframe later (via unnest) to put back into list
  return(as.numeric(ll))
}

# get_log_lik_wrapper <- function(b,a,nu,beta,phi,sigma,modelID='M0'){
#   th = c(b,a,nu,beta,phi,sigma)
#   producers = producers_list[[modelID]]
#   blocked = blocked_list[[modelID]]
#   ll = get_log_lik(overexpression_data,nTestOE,th,producers,blocked,times)
#   #easier for getting into right form in dataframe later (via unnest) to put back into list
#   return(lapply(1:nTestOE, FUN = function(i) ll[i,]))
# }

check_and_load_path <- function(fit_path){
  if (file.exists(fit_path)){
    estimates <- readRDS(fit_path)
  } else {
    #try on scratch
    scratch_path = paste('/scratch/harrison/FISH_data/old_fits/',stringr::str_split(fit_path,'/')[[1]][2],sep='')
    if (file.exists(scratch_path)){
      estimates <- readRDS(scratch_path)
    } else {
      stop('file does not exist, please run analysis or check a different computer')
    }
  }
  return(estimates)
}

# #informal test of log likelihood
# th = c(0.2, 10.0,  0.9, 1.0,  0.3)
# sig = rep(1.3,16)
# forward_simulate_wrapper(0.2,10,0.9,0.01,0.3,1.3,'M0')
# get_log_lik(overexpression_data,nTestOE,th,producers_list[['M0']],blocked_list[['M0']],times)
# alternative_log_lik_wrapper(10,0.2,0.9,1,0.3,1.3,'M0')


#get results from fitting simple model
cat(paste('\nloading fitted simple model parameters from model id: ',identifier_simple,'\n',sep=''))
path_simple = paste('fits/mrna_transport_estimates',identifier_simple,'.rds',sep='')
estimates_simple <- check_and_load_path(path_simple)
cat(paste('\nloading fitted density dependent model parameters from model id: ',identifier_DD,'\n',sep=''))
path_DD = paste('fits/mrna_transport_estimates',identifier_DD,'.rds',sep='')
estimates_DD <- check_and_load_path(path_DD)
#################

#get posterior parameters
draws_simple <- estimates_simple %>%
  spread_draws(a,b,nu,phi,sigma[i]) %>% 
  filter(.iteration<=num_draws) %>%
  spread(i,sigma,sep='sigma') %>%
  select(-.chain,-.iteration,-.draw) %>%
  mutate(beta = -1) 

draws_DD <- estimates_DD %>%
  spread_draws(a,b,nu,beta,phi,sigma[i]) %>% 
  filter(.iteration<=num_draws) %>%
  spread(i,sigma,sep='sigma') %>%
  select(-.chain,-.iteration,-.draw)
#################
# 
# log_lik_simple <- draws_simple %>%
#   mutate(modelID = list(models_to_sim[!(models_to_sim %in% models_with_density_dependence)])) %>%
#   unnest() %>%
#   mutate(time=list(times$ts4),
#          ll = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
#                                  get_log_lik_wrapper)) %>%
#   unnest(ll,time,.preserve=modelID) %>%
#   mutate(cellID = list(seq_len(16))) %>%
#   unnest(ll,cellID,.preserve=modelID) %>%
#   unnest()

###########

alternative_analysis <- function(mID,drop_egg_ind=FALSE,drop_cell_ind=FALSE){
  if (mID %in% models_with_density_dependence){
    q <- draws_DD
  } else {
    q <- draws_simple
  }
  q %<>% mutate(modelID = mID,
                drop_egg_ind=drop_egg_ind,drop_cell_ind=drop_cell_ind) %>%
    mutate(ll = purrr::pmap(list(b,a,nu,beta,phi,
                                 isigma1,isigma2,isigma3,isigma4,isigma5,isigma6,isigma7,isigma8,
                                 isigma9,isigma10,isigma11,isigma12,isigma13,isigma14,isigma15,isigma16,
                                 modelID,drop_egg_ind,drop_cell_ind),
                            alternative_log_lik_wrapper))
  alternative_log_lik_simple <- mapply(q[['ll']], FUN = function(x) as.numeric(x)) %>% t()
  loo1 <- loo(alternative_log_lik_simple)
  return(loo1)
}
get_model_weights <- function(models_to_sim,models_with_density_dependence,i,j){
  #i is the index of egg chamber to drop, j is index of 
aM_all <- list()
simple_models <- models_to_sim[!(models_to_sim %in% models_with_density_dependence)]
expose_stan_functions('mrna_transport_with_blocking.stan')
aM_all[simple_models] <- purrr::map(simple_models, function(mID) alternative_analysis(mID,drop_egg_ind = i,drop_cell_ind = j))
expose_stan_functions('mrna_transport_density_dependent_with_blocking.stan') #should replace previous functions and work with previous wrappers
DD_models <- models_to_sim[models_to_sim %in% models_with_density_dependence]
aM_all[DD_models] <- purrr::map(DD_models, function(mID) alternative_analysis(mID,drop_egg_ind = i,drop_cell_ind = j))
# purrr::map(aM_all,plot)

#get nicely formatted output from comparison
w_stack = loo::loo_model_weights(aM_all,method='stacking')
w_pseudo_bma = loo::loo_model_weights(aM_all,method="pseudobma")
return(list(w_stack = w_stack,w_pseudo_bma = w_pseudo_bma ))
}

z = get_model_weights(models_to_sim,models_with_density_dependence,0,0)
inhomogeneous_weights_df <- data_frame(model=models_to_sim,pseudo_bma=z$w_pseudo_bma[models_to_sim],stacking=z$w_stack[models_to_sim])
inhomogeneous_weights_df <- data_frame(model=names(z$w_pseudo_bma),pseudo_bma=as.numeric(z$w_pseudo_bma),stacking=as.numeric(z$w_stack))

method_names <- c(
  `pseudo_bma` = "Pseudo BMA+",
  `stacking` = "Stacking"
)
inhomogeneous_weights_df %>%
  gather(value=weight,key=method,-model) %>%
  ggplot(aes(x=model,y=weight)) +
  geom_col() +
  facet_wrap(~method, labeller = as_labeller(method_names)) + 
  theme_bw() + 
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 12)) +
  xlab('Model') + 
  ylab('Weight') +
ggsave(paste('plots/weights_model_comparison_',identifier_simple,'.eps',sep=''),device=cairo_ps)
