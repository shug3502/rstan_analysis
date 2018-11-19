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
  ###############
  th = c(b,a,nu,beta,phi,sigma)
  producers = producers_list[[modelID]]
  blocked = blocked_list[[modelID]]
  s = forward_simulate(nTestOE,th,producers,blocked,times)
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

#informal test of log likelihood
th = c(0.2, 10.0,  0.9, 1.0,  0.3,  1.3)
forward_simulate_wrapper(0.2,10,0.9,0.01,0.3,1.3,'M0')
get_log_lik(overexpression_data,nTestOE,th,producers_list[['M0']],blocked_list[['M0']],times)
alternative_log_lik_wrapper(10,0.2,0.9,1,0.3,1.3,'M0')


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
  spread_draws(a,b,nu,beta,phi,sigma) %>% 
  filter(.iteration<=num_draws) %>%
  select(b,a,nu,beta,phi,sigma)

#################

# post_pred_simple <- draws_simple %>%
#   mutate(modelID = list(models_to_sim[!(models_to_sim %in% models_with_density_dependence)])) %>%
#   unnest() %>%
#   mutate(time=list(times$ts4),
#          y_pred_oe = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
#                                  forward_simulate_wrapper)) %>%
#   unnest(y_pred_oe,time,.preserve=modelID) %>%
#   mutate(cellID = list(seq_len(16))) %>%
#   unnest(.preserve=modelID) %>%
#   unnest() %>%
#   group_by(time,cellID,modelID) %>%
#   mutate(y_lower=quantile(y_pred_oe,0.025),
#          y_upper=quantile(y_pred_oe,0.975),
#          y_median=quantile(y_pred_oe,0.5)) %>%
#   ungroup()
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

# out <- list()
# for (i in seq_len(nTestOE+1)){
#   out[[i]] <- list()
#   for (j in seq_len(17)){
#     print(c(i,j))
#     out[[i]][[j]] <- get_model_weights(models_to_sim,models_with_density_dependence,i-1,j-1)
#   }
# }
# print(out)
z = get_model_weights(models_to_sim,models_with_density_dependence,0,0)
weights_df <- data_frame(model=models_to_sim,pseudo_bma=z$w_pseudo_bma[models_to_sim],stacking=z$w_stack[models_to_sim])

method_names <- c(
  `pseudo_bma` = "Pseudo BMA+",
  `stacking` = "Stacking"
)
weights_df %>%
  gather(value=weight,key=method,-model) %>%
  ggplot(aes(x=model,y=weight)) +
  geom_col() +
  facet_wrap(~method, labeller = as_labeller(method_names)) + 
  theme_bw() + 
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) +
  xlab('Model') + 
  ylab('Weight') +
ggsave(paste('plots/weights_model_comparison_',identifier_simple,'.eps',sep=''),device=cairo_ps)

###################
# expose_stan_functions('mrna_transport_density_dependent_with_blocking.stan') #should replace previous functions and work with previous wrappers
# 
# post_pred_DD <- draws_DD %>%
#   mutate(modelID = list(models_to_sim[models_to_sim %in% models_with_density_dependence])) %>%
#   unnest() %>%
#   mutate(time=list(times$ts4),
#          y_pred_oe = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
#                                  forward_simulate_wrapper)) %>%
#   unnest(y_pred_oe,time,.preserve=modelID) %>%
#   mutate(cellID = list(seq_len(16))) %>%
#   unnest(.preserve=modelID) %>%
#   unnest() %>%
#   group_by(time,cellID,modelID) %>%
#   mutate(y_lower=quantile(y_pred_oe,0.025),
#          y_upper=quantile(y_pred_oe,0.975),
#          y_median=quantile(y_pred_oe,0.5)) %>%
#   ungroup()
# 
# log_lik_DD <- draws_DD %>%
#   mutate(modelID = list(models_to_sim[models_to_sim %in% models_with_density_dependence])) %>%
#   unnest() %>%
#   mutate(time=list(times$ts4),
#          ll = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
#                           get_log_lik_wrapper)) %>%
#   unnest(ll,time,.preserve=modelID) %>%
#   mutate(cellID = list(seq_len(16))) %>%
#   unnest(ll,cellID,.preserve=modelID) %>%
#   unnest()
# 
# #####################
# 
# full_log_lik_df <- full_join(log_lik_simple, log_lik_DD)
# 
# log_lik_df2 <- full_log_lik_df %>%
# #  filter(cellID!=1) %>%
#   group_by(modelID,b,a,nu,phi,sigma) %>%
#   summarise(ll_sum = sum(ll))
# 
# ######################
# 
# #add observed overexpression data
# xdata <- data.frame(rna = as.vector(overexpression_data),cellID = as.vector(matrix(rep(1:16,nTestOE),nrow=nTestOE,byrow=TRUE)),time = rep(times$ts4,16))
# observed_and_predictions <- full_join(full_join(post_pred_simple, post_pred_DD), xdata)
# 
# observed_and_predictions %>%
#   filter(round(time)==28) %>%
#   ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
#   geom_ribbon(alpha=0.4) +
#   geom_line() +
#   geom_point(aes(x=cellID,y=rna)) +
#   facet_wrap(~factor(modelID))
# 
# observed_and_predictions %>%
#   filter(modelID=='M6') %>%
#   ggplot(aes(x=cellID,y=y_median,ymin=y_lower,ymax=y_upper)) +
#   geom_ribbon(alpha=0.4) +
#   geom_line() +
#   geom_point(aes(x=cellID,y=rna)) +
#   facet_wrap(~factor(time),scales='free_y')  
# 
# g <- log_lik_df2 %>%
#   ggplot(aes(ll_sum,color=modelID)) +
#   geom_density() + 
#   scale_x_continuous(limits=c(-1500,-800))
#   print(g)
#   ggsave(paste('plots/compare_production_models_fwd_sim_',identifier_simple,'.eps',sep=''),device=cairo_ps)
#   
#   log_lik_df2 %>%
#     filter(as.numeric(substr(modelID,2,3))>10 | modelID=='M6' | modelID=='M3') %>%
#     ggplot(aes(ll_sum,color=modelID)) +
#     geom_density() + 
#     scale_x_continuous(limits=c(-1100,-850))
#   
#   ##########################
#   #plot log likelihood by cellID
#   h1 <- full_log_lik_df %>%
#     filter(cellID!=1) %>%
#     group_by(modelID,b,a,nu,phi,sigma,cellID) %>%
#     summarise(sum_ll = sum(ll)) %>%  #summing over egg chambers
#     group_by(modelID,cellID) %>% 
#     summarise(log_lik = median(sum_ll)) %>% #take median to summarise posterior
#     filter(modelID!='M10') %>% #model 10 is very bad
#     ggplot(aes(x=cellID,y=log_lik,color=modelID,group=modelID)) + 
#     geom_line() +
#     geom_point()
#   print(h1)
#   
#   #and by egg chamber
#   h2 <- full_log_lik_df %>%
#     filter(cellID!=1) %>%
#     group_by(modelID,b,a,nu,phi,sigma,time) %>%
#     summarise(sum_ll = sum(ll)) %>%  #summing over cells
#     group_by(modelID,time) %>% 
#     summarise(log_lik = median(sum_ll)) %>% #take median to summarise posterior
#     filter(modelID!='M10') %>% #model 10 is very bad
#     ggplot(aes(x=time,y=log_lik,color=modelID,group=modelID)) + 
#     geom_line() +
#     geom_point()
#   print(h2)
#   
#   h3 <- full_log_lik_df %>%
#     filter(cellID!=1) %>%
#     filter(time>5) %>%
#     group_by(modelID,b,a,nu,phi,sigma,cellID) %>%
#     summarise(sum_ll = sum(ll)) %>%  #summing over egg chambers
#     group_by(modelID,b,a,nu,phi,sigma) %>% 
#     filter(modelID!='M10') %>%
#     summarise(log_lik = sum(sum_ll)) %>% #sum again
#     group_by(modelID) %>%
#     summarise(mu = median(log_lik)) %>%
#     ggplot(aes(x=modelID,y=mu,color=modelID)) + 
#     geom_point()
#   print(h3)
#   
#   h4 <- full_log_lik_df %>%
#     group_by(modelID,b,a,nu,phi,sigma,cellID) %>%
#     summarise(sum_ll = sum(ll)) %>%  #summing over egg chambers
#     group_by(modelID,cellID) %>% 
#     filter(modelID!='M10') %>%
#     filter(as.numeric(substr(modelID,2,3))%%2==0) %>% 
#     ggplot(aes(x=factor(cellID),y=sum_ll,color=modelID)) + 
#     geom_violin()
#   print(h4)
#   
#   ###########################
#   # 
#   # analyse_given_model <- function(full_log_lik,nColumns,mID='M0'){
#   #   #mID is a string containing the model ID
#   #   log_lik_matrix = matrix(full_log_lik %>% filter(modelID==mID) %>% .$ll_sum,ncol=nColumns,byrow=FALSE)
#   #   loo1 <- loo(log_lik_matrix)
#   #   return(loo1) 
#   # }
# 
# #   #look at LOO by cellID
# #   log_lik_by_cellID <- full_log_lik_df %>%
# # #     filter(cellID!=1) %>% #if use this need 15 columns in matrix instead
# #     group_by(cellID,modelID,a,b,phi,sigma,nu,beta) %>%
# #     summarise(ll_sum=sum(ll))
# #   aM_cellID_list <- purrr::map(models_to_sim, function(x) analyse_given_model(log_lik_by_cellID,16,x))
# #   purrr::map(aM_cellID_list,plot)
# # 
# #   #look at LOO by egg chamber or time pt
# #   log_lik_by_time <- full_log_lik_df %>%
# #      filter(cellID!=1) %>%     
# #     group_by(time,modelID,a,b,phi,sigma,nu,beta) %>%
# #     summarise(ll_sum=sum(ll))
# #   aM_time_list <- purrr::map(models_to_sim, function(x) analyse_given_model(log_lik_by_time,nTestOE,x))
# #   purrr::map(aM_time_list,plot)
# #   
# #   ############################
# # #  try loo on all the data points
# #   log_lik_all <- full_log_lik_df %>%
# #     mutate(ll_sum=ll) 
# #   analyse_given_model(log_lik_all,16*nTestOE,'M0')
# #   
# # aM_all_list <- purrr::map(models_to_sim, function(x) analyse_given_model(log_lik_all,16*nTestOE,x))
# #   
#   