#setwd('~/Documents/FISH_data/rstan_analysis')
mrna_transport_inference_input_B <- function(identifier='v099',use_real_data=TRUE,run_mcmc=FALSE,nSamples=15,nTest=5,
                                      parametersToPlot = c("theta","phi","sigma"),verbose=FALSE,compare_via_loo=TRUE,
                                      show_diagnostic_plots=FALSE){
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# identifier = 'v053' #run identifier
# use_real_data <- TRUE
# run_mcmc <- TRUE
# nSamples = 15 #how many egg chambers segmented
# nTest = 5 #how many to test on

#############################################################
#Fit linear model to log(length) of egg chambers to get time of development.
#Could use area or volume of egg chamber. Compare to Jia 2016 Sci Reports.
#Area may be better as in some examples egg chambers are quite squashed.
#To get lengths read length.txt file in each folder

source('extract_times_and_scaling.R')
times = extract_times_and_scaling(nSamples,nTest)
#############################################################
m0 = c(0, rep(1,15)) #initial condition
th = c(6.8,132.8)
sig = 1.080
phi = 0.23

#############################################################
# source('get_nc_transition_matrix.R')
# B1 = get_nc_transition_matrix(0) %>% as.vector
# B2 = get_nc_transition_matrix(1) %>% as.vector

mc <- stan_model('model_comparison_at_stst_with_decay.stan')  #('mrna_transport5.stan')
expose_stan_functions(mc)
nu = 0.95
gamma = 0.01
B = construct_matrix(nu, gamma) %>% as.vector
#############################################################

if (use_real_data){
  #data taken from spot detection on 3 egg chambers
  if (verbose) print('using real data \n')
  #system('python ../custom_analysis/process_all_NCs.py',wait=TRUE)
  data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
  raw_data = data[times$sort_indices1,] #need to sort time series and correspondingly reorder rows
  exp_data = raw_data
  exp_data[is.na(exp_data)]=0 #stan can't deal with NAs
  #test_data = data[nSamples+sort_indices2,]
  test_data = data[times$sort_indices2,]
  } else {
  #sample from the model to get fake data
  print('using fake generated data')
  samples <- stan(file = 'mrna_transport6.stan',
                  data = list (
                    T  = nSamples+nTest,
                    y0 = m0,
                    t0 = times$t0,
                    ts = times$ts2,
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
  
  s <- rstan::extract(samples,permuted=FALSE)
  #plot(s[1,1,seq(from=1, to=(nSamples*16-1), by=nSamples)])
  #plot(s[1,1,seq(from=nSamples, to=(nSamples*16), by=nSamples)])
  boxplot(s[,1,seq(from=nSamples, to=((nSamples+nTest)*16), by=nSamples)])
  test_data = matrix(s[1,1,1:(16*(nSamples+nTest))],nrow=(nSamples+nTest),byrow=FALSE) #this is our fake data
  exp_data = test_data[!(times$ts2 %in% times$ts3),] # for plotting
}


##########################
if (run_mcmc) {
estimates <- stan(file = 'mrna_transport5.stan',
                  data = list (
                    y = exp_data,
                    T1  = nSamples,
                    T2 = nSamples+nTest,
                    y0 = m0,
                    t0 = times$t0,
                    ts1 = times$ts1,
                    ts2 = times$ts2,
                    nu = c(nu)
                  ),
                  seed = 42,
                  chains = 4,
		warmup = 500,
		iter = 1000
)

tryCatch({
  #estimates@stanmodel@dso <- new("cxxdso") #seems have to do this to be able to save :S
  saveRDS(estimates, file = paste('fits/mrna_transport_estimates',identifier,'.rds',sep='')) # to load use readRDS("fit.rds")
}, warning = function(war) {  
  # warning handler picks up where error was generated
  print(paste("WARNING:  ",war))
}, error = function(err) {  
  # error handler picks up where error was generated
  print(paste("NOT SAVED ERROR:  ",err))
}) # END tryCatch
} else {
  estimates = readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep=''))
}

#parametersToPlot = c("theta","phi","sigma") #c('mu','tau','psi','zeta')# c('mu','tau','phi')
if (verbose) print(estimates, pars = parametersToPlot)

if (compare_via_loo){
  library(loo)
  log_lik_1 <- extract_log_lik(estimates)
  loo_1 <- loo(log_lik_1)
  print(loo_1)
}  

#######################
#visualisation
###############################
#cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
pairs(estimates, pars = parametersToPlot)
#dev.off()

#look at posterior predictive distn
source('post_pred_plot.R')
p1 <- post_pred_plot(test_data,times$ts2,nTest+nSamples,'y_pred',estimates,identifier,title_stem='plots/posterior_pred',ts_test=times$ts3)

if (show_diagnostic_plots) {
  source('mcmcDensity.R')
  mcmcDensity(estimates, parametersToPlot, byChain = TRUE)
  ggsave(paste('plots/denisty',identifier, '.eps',sep=''),device=cairo_ps)
  
  library(bayesplot)
  draws <- as.array(estimates, pars=parametersToPlot)
  mcmc_trace(draws)
  ggsave(paste('plots/trace',identifier, '.eps',sep=''),device=cairo_ps)
  mcmc_intervals(draws,pars=c('a','b','phi'))
  ggsave(paste('plots/intervals',identifier, '.eps',sep=''),device=cairo_ps)
  color_scheme_set("purple")
  mcmc_areas(draws,pars=c('a','b','phi'))
  ggsave(paste('plots/areas',identifier, '.eps',sep=''),device=cairo_ps)
  color_scheme_set("brightblue")
  mcmc_scatter(draws,pars=c('a','b'))
  ggsave(paste('plots/scatter',identifier, '.eps',sep=''),device=cairo_ps)
  
  cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
  pairs(estimates, pars = parametersToPlot)
  dev.off()
}
######################
#cos saving wasn't working
e <- rstan::extract(estimates,pars=parametersToPlot,permuted=TRUE)
save(e,file=paste('fits/alt_save',identifier,'.Rsave',sep=''))

if (compare_via_loo){
  return(list(p1,loo_1,estimates))
} else {
  return(p1)
}
}
