#setwd('~/Documents/FISH_data/rstan_analysis')
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
identifier = 'MCv056' #run identifier
use_real_data <- TRUE
run_mcmc <- TRUE
nSamples = 15 #how many egg chambers segmented
nTest = 5 #how many to test on
nTotal = nSamples + nTest

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
sig = 150 
phi = 0.289 #careful can cause issues if phi is not 1 when looking at stst

#############################################################
my_normaliser <- function(dataset) apply(dataset,1,function(x)x/x[1]) %>% t()

if (use_real_data){
  #data taken from spot detection on 3 egg chambers
  print('using real data \n')
  #system('python ../custom_analysis/process_all_NCs.py',wait=TRUE)
  data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
  raw_data = data[times$sort_indices1,] #need to sort time series and correspondingly reorder rows
  exp_data = raw_data
  exp_data[is.na(exp_data)]=0 #stan can't deal with NAs
  full_data = data[times$sort_indices2,]
  normalised_data = exp_data %>% my_normaliser #divide by amount in oocyte to normalise for stst
  test_data = data[!(times$sort_indices2 %in% times$sort_indices1),]
  #test_data = apply(test_data,1,function(x)x/x[1]) %>% t()
} else {
  #sample from the model to get fake data
  print('using fake generated data')
  #############################################################
  source('get_nc_transition_matrix.R')
  B1 = get_nc_transition_matrix(0) %>% as.vector
  B2 = get_nc_transition_matrix(1) %>% as.vector
  
  mc <- stan_model('model_comparison_at_stst2.stan')
  expose_stan_functions(mc)
  nu = 0.72
  # set_oocyte_absorbing = function(B){
  #   B[,1] = 0
  #   return(B)
  # }
  #currently will not be using the absorbing oocyte fn
  B = construct_matrix(nu) %>% as.vector
  #############################################################
  samples <- stan(file = 'mrna_transport6.stan',
                  data = list (
                    T  = nTotal,
                    y0 = m0,
                    t0 = (times$t0)*10,
                    ts = (times$ts2)*10,
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
  plot(s[1,1,seq(from=1, to=(nTotal*16-1), by=nTotal)])
  plot(s[1,1,seq(from=nSamples, to=(nTotal*16), by=nTotal)])
  boxplot(s[,1,seq(from=nSamples, to=(nTotal*16), by=nTotal)])
  #   exp_data = matrix(s[1,1,1:(16*nSamples)],nrow=nSamples,byrow=FALSE) #this is our fake data
  full_data = matrix(s[1,1,1:(16*(nTotal))],nrow=(nTotal),byrow=FALSE) #this is our fake data
  exp_data = full_data[!(times$ts2 %in% times$ts3),] 
  normalised_data = exp_data %>% my_normaliser #divide by amount in oocyte to normalise for stst
  boxplot(normalised_data)
  test_data = full_data[(times$ts2 %in% times$ts3),]
  #test_data = apply(test_data,1,function(x)x/x[1]) %>% t()
}

##########################
if (run_mcmc) {
  estimates <- stan(file = 'model_comparison_at_stst2.stan',
                    data = list (
                      T1 = nSamples,
                      T2 = nTest,
                      y_obs = normalised_data
                    ),
                    seed = 42,
                    chains = 4,
                    warmup = 1000,
                    iter = 2000
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

parametersToPlot = c('nu','sigma','phi') #c('mu','tau','zeta','xi') 
print(estimates, pars = parametersToPlot)

source('hist_treedepth.R')
hist_treedepth(estimates)

#######################
#visualisation
###############################
cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
pairs(estimates, pars = parametersToPlot)
dev.off()
#ggsave(paste('plots/pairs',identifier, '.eps',sep=''),device=cairo_ps)

#look at posterior predictive distn
source('post_pred_plot_at_stst.R')
#post_pred_plot_at_stst((test_data %>% my_normaliser),times$ts3,nTest,'y_sim',estimates,identifier,title_stem='plots/posterior_pred_stst')
post_pred_plot_at_stst((full_data %>% my_normaliser),times$ts2,nTotal,'y_sim',estimates,identifier,title_stem='plots/posterior_pred_stst')
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
# mcmc_scatter(draws,pars=parametersToPlot)
# ggsave(paste('plots/scatter',identifier, '.eps',sep=''),device=cairo_ps)


######################
#cos saving wasn't working
e <- rstan::extract(estimates,pars=parametersToPlot,permuted=TRUE)
save(e,file=paste('fits/alt_save_MC',identifier,'.Rsave',sep=''))
