synthetic_mixture_model_test <- function(identifier = 'syntheticv001',
                                         run_mcmc=TRUE,
                                         parametersToPlot = c('nu','phi','xi'),
                                         nSamples = 3,
                                         nTest=1,
                                         nTestOE=0,
                                         rc_index_real = 16)
{
  ## trying mixture model with synthetic data to check for identifiability
  library(rstan)
  library(mvtnorm)
  library(dplyr)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  source('extract_times_and_scaling.R')
  times = extract_times_and_scaling(nSamples,nTest,nTestOE)
  #############################################################
  
  expose_stan_functions('model_comparison_at_stst3.stan')
  #suppose we observe data from one of the mixture components and feed this as data to the mixture model
  #can we identify the correct mixture component?
  nu_real = 0.95
  phi_real = 0.35
  xi_real = 0.1
  
  
  raw_data = get_k2(nu_real,get_RC_from_dict(rc_index_real))
  data = (raw_data * c(1,rep(1/phi_real,15))) %>% rep(.,nSamples) %>% matrix(.,ncol=16,byrow = TRUE) + matrix(abs(xi_real*rnorm(16*nSamples)),ncol=16)
  plot(data[1,])
  ##############################################################
  
  
  if (run_mcmc) {
    stan_file <- 'model_comparison_at_stst3.stan'
    estimates <- stan(file = stan_file,
                      data = list (
                        T1 = nSamples,
                        T2 = nTest+nTestOE,
                        y_obs = data
                      ),
                      seed = 42,
                      chains = 4,
                      warmup = 1000,
                      iter = 2000,
                      init = function() list(nu=nu_real, xi=xi_real, phi=phi_real),
                      control=list(adapt_delta=0.99)
    )
    
    tryCatch({
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
  print(estimates, pars = parametersToPlot)
  pairs(estimates, pars = parametersToPlot)
  
  source('post_pred_plot_at_stst.R')
  data <- rbind(data,(raw_data * c(1,rep(1/phi_real,15)))) #for single test examples
  p1 <- post_pred_plot_at_stst(data,times$ts2,nSamples+nTest+nTestOE,
                               'y_sim',estimates,identifier,title_stem='plots/posterior_pred_stst',
                               ts_test=vector(),OE_test=times$ts4) #times$ts3
  return(list(p1,estimates))
}

for (i in seq_len(1)){
  identifier = paste('syntheticv10',i,sep='')
  p <- synthetic_mixture_model_test(identifier, run_mcmc=FALSE, nSamples=10,  rc_index_real=i)
}
