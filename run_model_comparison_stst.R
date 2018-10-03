#setwd('~/Documents/FISH_data/rstan_analysis')
#write as function, may need to think about what to do if want to run from command line again
run_model_comparison_stst <- function(identifier='MCv099',use_real_data=TRUE,run_mcmc=FALSE,nSamples=9,nTest=11,nTestOE=9,
                                      parametersToPlot = c('nu','xi','phi'),verbose=FALSE,
                                      compare_via_loo=FALSE, show_diagnostic_plots=FALSE,
                                      use_mixture_model=FALSE, train_on_OE=FALSE){
  #steady state analysis only for the simplest model without decay or other complications
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('extract_times_and_scaling.R')
times = extract_times_and_scaling(nSamples,nTest,nTestOE)
#############################################################
m0 = c(0, rep(1,15)) #initial condition
th = c(6.8,132.8)
sig = 150 
phi = 0.289 #careful can cause issues if phi is not 1 when looking at stst

#############################################################
my_normaliser <- function(dataset) apply(dataset,1,function(x)x/x[1]) %>% t()
#########

if (use_real_data){
  #data taken from spot detection on 3 egg chambers
  if (verbose) print('using real data \n')
  #system('python ../custom_analysis/process_all_NCs.py',wait=TRUE)
  data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
  raw_data = data[times$sort_indices1,] #need to sort time series and correspondingly reorder rows
  exp_data = raw_data
  exp_data[is.na(exp_data)]=0 #stan can't deal with NAs
  overexpression_data = matrix(as.numeric(read.csv('data/exp_data_overexpression.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
  overexpression_data = overexpression_data[times$sort_indices4,] #need to sort time series and correspondingly reorder rows
  test_data = rbind(data,overexpression_data)
  test_data = test_data[times$sort_indices2,]
  test_data[is.na(test_data)]=0
  WT_test_data = data[times$sort_indices3,]
  WT_test_data[is.na(WT_test_data)]=0
  if (train_on_OE){
    swap <- function(name1, name2, envir = parent.env(environment()))
    {
      temp <- get(name1, pos = envir)
      assign(name1, get(name2, pos = envir), pos = envir)
      assign(name2, temp, pos = envir)
    }
    swap('exp_data','overexpression_data')
    swap('nSamples','nTestOE')
    temp <- list(times$ts1,times$ts4)
    times$ts1 <- temp[[2]]
    times$ts4 <- temp[[1]]
    print('swapped')
    print(times$ts1)
    print(nSamples)
  }
  normalised_data = exp_data %>% my_normaliser #divide by amount in oocyte to normalise for stst
  #for data driven prior from eyeballing which mixture egg chambers are from
  alpha = matrix(rep(1,(nSamples+nTestOE+nTest)*16),ncol=16)
  alpha[1,1] = 15
  alpha[2,16] = 15
  alpha[3,c(6,11,16)] = 13
  alpha[4,c(12,16)] = 14
  alpha[5,10] = 15
  alpha[6,c(10,11)] = 14
  alpha[7,c(1,16)] = 14
  alpha[8,c(6,9)] = 14
  alpha[9,c(6,11)] = 14 
} else {
  warning('TODO: update simulated data for full model')
  return(0)
}
############################
if (run_mcmc) {
  stan_file <- case_when(use_mixture_model ~ 'model_comparison_at_stst4.stan',
                        TRUE ~ 'model_comparison_at_stst2.stan')
  estimates <- stan(file = stan_file,
                    data = list (
                      T1 = nSamples,
                      T2 = nTest+nTestOE,
                      y_obs = normalised_data,
                      alpha = alpha
                    ),
                    seed = 42,
                    chains = 4,
                    warmup = 1000,
                    iter = 2000,
                    init = function() list(nu=0.9, xi=0.4, phi=0.3)
#                    control = list(adapt_delta = 0.99)
  )
  
  tryCatch({
    #estimates@stanmodel@dso <- new("cxxdso") #seems have to do this to be able to save :S
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

if (verbose) print(estimates, pars = parametersToPlot)

#######################
#visualisation
###############################
#cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
pairs(estimates, pars = parametersToPlot)
#dev.off()
#ggsave(paste('plots/pairs',identifier, '.eps',sep=''),device=cairo_ps)

#look at posterior predictive distn
source('post_pred_plot_at_stst.R')
#post_pred_plot_at_stst((test_data %>% my_normaliser),times$ts3,nTest,'y_sim',estimates,identifier,title_stem='plots/posterior_pred_stst')
p1 <- post_pred_plot_at_stst((test_data %>% my_normaliser),times$ts2,nSamples+nTest+nTestOE,
                             'y_sim',estimates,identifier,title_stem='plots/posterior_pred_stst',
                             ts_test=vector(),OE_test=times$ts4) #times$ts3

if (compare_via_loo){
  library(loo)
  log_lik_1 <- extract_log_lik(estimates)
  loo_1 <- loo(log_lik_1)
  print(loo_1)
} else {
  loo_1 <- NA
} 

if (show_diagnostic_plots) {
  source('hist_treedepth.R')
  hist_treedepth(estimates)
  source('mcmcDensity.R')
  mcmcDensity(estimates, parametersToPlot, byChain = TRUE)
  ggsave(paste('plots/denisty',identifier, '.eps',sep=''),device=cairo_ps)
  
  library(bayesplot)
  draws <- as.array(estimates, pars=parametersToPlot)
  mcmc_trace(draws)
  ggsave(paste('plots/trace',identifier, '.eps',sep=''),device=cairo_ps)
  mcmc_intervals(draws,pars=parametersToPlot)
  ggsave(paste('plots/intervals',identifier, '.eps',sep=''),device=cairo_ps)
  color_scheme_set("purple")
  mcmc_areas(draws,pars=parametersToPlot)
  ggsave(paste('plots/areas',identifier, '.eps',sep=''),device=cairo_ps)
  color_scheme_set("brightblue")
  cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
  pairs(estimates, pars = parametersToPlot)
  dev.off()
  # mcmc_scatter(draws,pars=parametersToPlot)
  # ggsave(paste('plots/scatter',identifier, '.eps',sep=''),device=cairo_ps)
}


######################
#cos saving wasn't working
e <- rstan::extract(estimates,pars=parametersToPlot,permuted=TRUE)
save(e,file=paste('fits/alt_save_MC',identifier,'.Rsave',sep=''))

return(list(p1,estimates,loo_1))
}

