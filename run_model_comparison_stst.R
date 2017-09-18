setwd('~/Documents/FISH_data/rstan_analysis')
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
identifier = 'MCv041' #run identifier
use_real_data <- FALSE
run_mcmc <- TRUE
nSamples = 5 #how many egg chambers segmented
nTest = 5 #how many to test on
nTotal = nSamples + nTest

#############################################################
#Fit linear model to log(length) of egg chambers to get time of development.
#Could use area or volume of egg chamber. Compare to Jia 2016 Sci Reports.
#Area may be better as in some examples egg chambers are quite squashed.
#To get lengths read length.txt file in each folder

# egg_chamber_areas <- rep(0,nTotal)
# #stages <- rep(0,nSamples) #don't need to extract estimated stages for each egg chamber example
# for (j in 1:(nTotal)){
#   egg_chamber_areas[j] <- as.numeric(read.table(paste('../data/Example',j,'/area.txt',sep='')))
# }
# 
# t0 = log(400)
# log_areas = sort.int(log(egg_chamber_areas[1:nSamples]),index.return=TRUE)
# ts1 = log_areas$x #ts needs to be an ordered time vector
# log_areas_test = sort.int(log(egg_chamber_areas),index.return=TRUE)
# ts2 = log_areas_test$x #includes the extra test sets or ts = setdiff(ts2,ts1)
# sort_indices1 = log_areas$ix
# sort_indices2 = log_areas_test$ix
# ts3 = setdiff(ts2,ts1) 
# log_areas3 = sort.int(log(egg_chamber_areas[(nSamples+1):(nSamples+nTest)]),index.return=TRUE)
# sort_indices3 = log_areas3$ix
source('extract_times_and_scaling.R')
times = extract_times_and_scaling(nSamples,nTest)
#############################################################
m0 = c(0, rep(1,15)) #initial condition
th = c(6.8,132.8)
sig = 0.10 
phi = 0.289

#############################################################
source('get_nc_transition_matrix.R')
B1 = get_nc_transition_matrix(0) %>% as.vector
B2 = get_nc_transition_matrix(1) %>% as.vector

mc <- stan_model('model_comparison_test.stan')
expose_stan_functions(mc)
nu = 0.90
set_oocyte_absorbing = function(B){
  B[,1] = 0
  return(B)
}
#currently will not be using the absorbing oocyte fn
B = construct_matrix(nu) %>% as.vector
#############################################################

if (use_real_data){
  #data taken from spot detection on 3 egg chambers
  print('using real data \n')
  #system('python ../custom_analysis/process_all_NCs.py',wait=TRUE)
  data = matrix(as.numeric(read.csv('../data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
  raw_data = data[times$sort_indices1,] #need to sort time series and correspondingly reorder rows
  exp_data = raw_data
  exp_data[is.na(exp_data)]=0 #stan can't deal with NAs
  test_data = data[times$sort_indices2,]
} else {
  #sample from the model to get fake data
  print('using fake generated data')
  samples <- stan(file = 'mrna_transport6.stan',
                  data = list (
                    T  = nTotal,
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
  plot(s[1,1,seq(from=1, to=(nTotal*16-1), by=nTotal)])
  plot(s[1,1,seq(from=nSamples, to=(nTotal*16), by=nTotal)])
  boxplot(s[,1,seq(from=nSamples, to=(nTotal*16), by=nTotal)])
  #   exp_data = matrix(s[1,1,1:(16*nSamples)],nrow=nSamples,byrow=FALSE) #this is our fake data
  test_data = matrix(s[1,1,1:(16*(nTotal))],nrow=(nTotal),byrow=FALSE) #this is our fake data
  exp_data = test_data[!(times$ts2 %in% times$ts3),] 
  normalised_data = apply(exp_data,1,function(x)x/x[1]) %>% t() #divide by amount in oocyte to normalise for stst
}


##########################
if (run_mcmc) {
  estimates <- stan(file = 'model_comparison_at_stst.stan',
                    data = list (
                      T1 = nSamples,
                      T2 = nTest,
                      y_obs = normalised_data
                    ),
                    seed = 42,
                    chains = 4,
                    warmup = 500,
                    iter = 1000,
                    control = list(adapt_delta = 0.9, 
                                   max_treedepth = 15)
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

parametersToPlot = c('mu','tau','zeta','xi')
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
post_pred_plot_at_stst(exp_data,times$ts1,nSamples,'y_sim',estimates,identifier,title_stem='plots/posterior_pred_stst')
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
mcmc_scatter(draws,pars=parametersToPlot)
ggsave(paste('plots/scatter',identifier, '.eps',sep=''),device=cairo_ps)


######################
#cos saving wasn't working
e <- rstan::extract(estimates,pars=parametersToPlot,permuted=TRUE)
save(e,file=paste('fits/alt_save_MC',identifier,'.Rsave',sep=''))
