#Aims to get predictions and data on RNA accumulation in the oocyte
#From forward simulate.
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
identifier_simple = 'v470_with_outliers_M0_simple' #'v431WT_simple' #load fitted posterior parameters from this model
identifier_DD = 'v470_with_outliers_M8_density_dependent_4a'
nSamples = 20
nTestOE = 14
num_draws = 100
gamma = 2
models_to_sim = paste('M',seq(from=0,to=7),sep='')
models_with_blocking = c('M1','M5','M6','M7')
models_with_density_dependence = c('M2','M4','M6','M7')
models_with_data_driven_production = c('M3','M4','M5','M7')
times = extract_times_and_scaling(nSamples,0,nTestOE)
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

# #######################

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

#get results from fitting simple model
cat(paste('\nloading fitted simple model parameters from model id: ',identifier_simple,'\n',sep=''))
path_simple = paste('fits/mrna_transport_estimates',identifier_simple,'.rds',sep='')
estimates_simple <- check_and_load_path(path_simple)
cat(paste('\nloading fitted density dependent model parameters from model id: ',identifier_DD,'\n',sep=''))
path_DD = paste('fits/mrna_transport_estimates',identifier_DD,'.rds',sep='')
estimates_DD <- check_and_load_path(path_DD)
# #################
# 
#get posterior parameters
draws_simple <- estimates_simple %>%
  spread_draws(a,b,nu,phi,sigma) %>%
  filter(.iteration<=num_draws) %>%
  select(b,a,nu,phi,sigma) %>%
  mutate(beta=-1)
# 
post_pred_simple <- draws_simple %>%
  mutate(modelID = list(models_to_sim[!(models_to_sim %in% models_with_density_dependence)])) %>%
  unnest() %>%
  mutate(phenotype=list(c("WT","OE"))) %>%
  unnest() %>%
  mutate(a = ifelse(phenotype=="OE",gamma*a,a),time=ifelse(phenotype=="OE",list(times$ts4),list(times$ts1))) %>%
  mutate(y_pred_oe = purrr::pmap(list(b,a,nu,beta,phi,sigma,modelID),
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
# 
# ###############
# 
# expose_stan_functions('mrna_transport_density_dependent_with_blocking.stan') #should replace previous functions and work with previous wrappers
# 
# draws_DD <- estimates_DD %>%
#   spread_draws(a,b,nu,beta,phi,sigma) %>% 
#   filter(.iteration<=num_draws) %>%
#   select(b,a,nu,beta,phi,sigma)
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
# #add observed overexpression data
# xdata <- data.frame(rna = as.vector(overexpression_data),cellID = as.vector(matrix(rep(1:16,nTestOE),nrow=nTestOE,byrow=TRUE)),time = rep(times$ts4,16))
# observed_and_predictions <- full_join(full_join(post_pred_simple, post_pred_DD), xdata)
# 
# fig_df <- observed_and_predictions %>%
#   filter(cellID==1) %>%
#   group_by(time,modelID) %>%
#   summarise(number = unique(rna),
#             model_number = unique(y_median))
# 
# 
# h1 <- ggplot(data = fig_df %>% filter(modelID=='M0'), aes(x=time,y=number)) + 
#   geom_point() +
#   geom_smooth(method = lm, se = TRUE) +
#   theme_bw() +
#   theme(text = element_text(size = 12), axis.text = element_text(size = 12),
#         legend.position = "none", strip.text = element_text(size = 8)) +
#   xlab('Time') + 
#   ylab('mRNA complexes in oocyte')
# print(h1)
# 
# h2 <- ggplot(data = fig_df %>% filter(modelID=='M0'), aes(x=time,y=model_number,color=modelID)) + 
#   geom_point() +
#   geom_smooth(method = lm, se = TRUE) +
#   theme_bw() +
#   theme(text = element_text(size = 12), axis.text = element_text(size = 12),
#         strip.text = element_text(size = 8)) +
#   xlab('Time') + 
#   ylab('mRNA complexes in oocyte')
# print(h2)
# 
# library(patchwork)
# h1+h2
# ################
# fig_df %>%
#   gather(sim, rna,-time,-modelID) %>%
#   ggplot(aes(x=time,y=rna,color=modelID)) + 
#   geom_point() +
#   geom_smooth(method=lm) +
#   facet_wrap(~sim) + 
#   theme_bw()
# 
# fig_df %>%
#   gather(sim, rna,-time,-modelID) %>%
#   filter(modelID=='M0') %>%
#   ggplot(aes(x=time,y=rna,color=sim)) + 
#   geom_point() +
#   geom_smooth(method=lm) +
#   theme_bw()
# 
# fig_df %>%
#   gather(sim, rna,-time,-modelID) %>%
#   group_by(modelID) %>%
#   mutate(time=15*time) %>%
# do(tidy(lm(rna ~ time, data = .)))
# ##really confused by the answers from this, slightly less than 4. Similar to the optimal values for gamma
# ##But this should give estimates of gamma*a
# 
# observed_and_predictions %>%
#   filter(cellID==1) %>%
#   group_by(time,modelID) %>%
#   summarise(number = unique(rna),
#             model_number = unique(y_median),
#             temp = median(temp)) %>%
#   gather(sim, rna,-time,-modelID) %>%
#   filter(modelID=='M0') %>%
#   ggplot(aes(x=time,y=rna,color=sim)) + 
#   geom_point() +
#   geom_smooth(method=lm) +
#   theme_bw()
# ##################
# 
# y_analytic_pred <- draws_simple %>%
#   mutate(time = list(times$ts4)) %>%
#   unnest() %>%
#   mutate(y_anal = purrr::pmap(list(b,2*a,nu,time),analytic_soln),
#          cellID=list(seq(16))) %>%
#   unnest()
# y_df <- y_analytic_pred %>%
#   group_by(time,cellID) %>%
#   summarise(y_med = median(y_anal)*0.345)
# 
# ggplot(y_df %>% filter(cellID==1),aes(x=time,y=y_med)) +
#   geom_point() +
#   geom_smooth(method=lm)
# 
# y_df %>% filter(cellID==1) %>% ungroup() %>% mutate(time=15*time) %>% do(tidy(lm(y_med ~ time, data = .)))

###############
#with analytic model predictions and data on same plot
source('analytic_soln.R')
#gamma should be multiplier for overexpresser production
y_analytic_wt <- draws_simple %>%
  mutate(phenotype='WT',time = list(times$ts1)) %>%
  unnest() %>%
  mutate(y_analytic = purrr::pmap(list(b,a,nu,time),analytic_soln),
         cellID=list(seq(16))) %>%
  unnest()
y_analytic_oe <- draws_simple %>%
  mutate(phenotype='OE',time = list(times$ts4)) %>%
  unnest() %>%
  mutate(y_analytic = purrr::pmap(list(b,gamma*a,nu,time),analytic_soln),
         cellID=list(seq(16))) %>%
  unnest()
y_df <- full_join(y_analytic_oe,y_analytic_wt) %>%
  group_by(time,cellID,phenotype) %>%
  summarise(Model = median(y_analytic*if_else(cellID==1,phi,1)),
            model_lower=quantile(y_analytic*if_else(cellID==1,phi,1),0.025),
            model_upper=quantile(y_analytic*if_else(cellID==1,phi,1),0.975))

######## 
#add data
#add observed overexpression data
data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
wildtype_data = data[times$sort_indices1,] #need to sort time series and correspondingly reorder rows
oe_data <- data.frame(Data = as.vector(overexpression_data),cellID = as.vector(matrix(rep(1:16,nTestOE),nrow=nTestOE,byrow=TRUE)),time = rep(times$ts4,16),phenotype='OE')
wt_data <- data.frame(Data = as.vector(wildtype_data),cellID = as.vector(matrix(rep(1:16,nSamples),nrow=nSamples,byrow=TRUE)),time = rep(times$ts1,16),phenotype='WT')
xdata = full_join(oe_data,wt_data)
observed_and_predictions <- full_join(y_df, xdata) %>% gather(method,rna,-time,-cellID,-phenotype,-model_lower,-model_upper)

###########
#plot
ggplot(observed_and_predictions, aes(x=time,y=rna,color=method,shape=phenotype)) +
  geom_point() + 
  facet_wrap(~cellID) +
  theme_bw()

#view as difference between WT and OE
ggplot(observed_and_predictions %>% filter(cellID==1), aes(x=time,y=rna)) +
  geom_point(aes(color=phenotype,shape=phenotype)) +
  geom_ribbon(aes(ymin=model_lower,ymax=model_upper,color=phenotype),alpha=0.3) +
  facet_wrap(~method) +
  theme_bw()

df <- post_pred_simple %>%
  group_by(time,cellID,phenotype,modelID) %>%
  summarise(rna=unique(y_median),method=unique(modelID)) %>%
  full_join(observed_and_predictions)

ggplot(df %>% filter(cellID==1), aes(x=time,y=rna,color=phenotype, shape=method)) +
  geom_point() + 
  theme_bw()

############
bl = sapply(1:14, function(n) paste(rep(" ",n),collapse=""))
egg_chamber_layout_levels <- c(bl[1:2], 4, bl[3:4],
                               bl[5], 8, 6, 2, bl[6],
                               bl[7], 12, 7, 3, bl[8],
                               16, 14, 10, 5, 1,
                               bl[9], 15, 11, 9, bl[10],
                               bl[11:12], 13, bl[13:14])
observed_and_predictions$cellID <- factor(observed_and_predictions$cellID,
                                          levels = egg_chamber_layout_levels) 

h <- ggplot( observed_and_predictions %>% filter(phenotype=='WT'), aes(x=time,y=rna)) +
  geom_point(aes(color=method,shape=method)) +
  geom_ribbon(aes(ymin=model_lower,ymax=model_upper,color=method),alpha=0.3) +
  facet_wrap(~cellID, ncol = 5, drop = F, strip.position="bottom") +
  theme_classic() + 
  theme(strip.background=element_blank())

