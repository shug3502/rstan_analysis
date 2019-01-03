#Fits three types of model, indicated by the model_str argument:
# simple - basic model given in the manuscript without decay
# decay - simple model with decay
# density_dependent - density dependent model with hill type functional form with n=1
#For comparison with OE data, use the forward simulations based on parameter values from fitting these models to WT data.
#See forward_simulate.R for this

#JH last updated 06/11/18 

##################

mrna_transport_inference_full <- function(identifier='full_v099',use_real_data=TRUE,run_mcmc=FALSE,nSamples=15,nTest=5,nTestOE=3,
                                          parametersToPlot = c("theta","phi","sigma","a","b"),verbose=FALSE,compare_via_loo=FALSE,
                                          show_diagnostic_plots=FALSE, train_on_OE=FALSE, alpha=2, model_str='simple'
                                          ){
  library(rstan)
  library(mvtnorm)
  library(dplyr)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

stopifnot(any(model_str %in% c('simple','decay','density_dependent'))) #check if model is one of simple, decay or DD
#check folder structure
stopifnot(dir.exists('fits') & dir.exists('plots'))
  
  #############################################################
  #Fit linear model to log(length) of egg chambers to get time of development
  #using area of egg chamber. Compare to Jia 2016 Sci Reports.
  
  source('get_just_areas.R')
  egg_chamber_areas <- get_just_areas(nSamples,nTest,nTestOE)
  stages <- get_just_stages(nSamples,nTest,nTestOE)
  time_hrs <- convert_to_hrs(stages)
  source('extract_times_and_scaling.R')
  times = extract_times_and_scaling(nSamples,nTest,nTestOE)
  m0 = c(0, rep(0,15)) #initial condition
  #############################################################
  
  #############################################################
  
  if (use_real_data){
    #data taken from spot detection on egg chambers via FISH quant
    if (verbose) print('using real data \n')
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
      if (verbose){
        print('swapped')
        print(times$ts1)
        print(nSamples)
      }
    }
  } else {
    warning('TODO: update simulated data for full model')
    #sample from the model to get fake data 
    mc <- stan_model('mrna_transport_with_blocking.stan') 
    expose_stan_functions(mc)
    th = c(0.2,10)
    sig = 1.25
    phi = 0.345
    nu = 0.92
    gamma = 0.3
    B = construct_matrix(nu, gamma) %>% as.vector

    print('using fake generated data')
    samples <- stan(file = 'mrna_transport6.stan',
                    data = list (
                      T  = nSamples+nTest+nTestOE,
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
    test_data = matrix(s[1,1,1:(16*(nSamples+nTest+nTestOE))],nrow=(nSamples+nTest+nTestOE),byrow=FALSE) #this is our fake data
    exp_data = test_data[!(times$ts2 %in% times$ts3),] 
    overexpression_data = test_data[(times$ts2 %in% times$ts4),]
  }

    #alpha is multiplier for how much more rna is produced in the overexpression mutant. Estimate via 'how_much_more_produced_in_oe.stan'
    producers = matrix(rep(alpha,nTestOE*16),ncol=16)
    producers[,1] = 0 #homogeneous production
    #data driven producers are considered in the forward_simulate.R
    blocked = matrix(rep(1,nTestOE*3),ncol=3) #this will make sure nothing is blocked in posterior predictive sims
    #similarly blocking is considered there too

  ##########################
  if (run_mcmc) {
    # if (use_prior_predictive){
    # stan_file = case_when( model_str == 'simple' ~ 'prior_predictive_with_blocking.stan',
    #                        model_str == 'decay' ~ 'prior_predictive_with_decay.stan',
    #                        model_str == 'density_dependent' ~ 'prior_predictive_density_dependent_with_blocking.stan'
    # )
    # } else {
     stan_file = case_when( model_str == 'simple' ~ 'mrna_transport_with_blocking.stan',
                             model_str == 'decay' ~ 'mrna_transport_with_decay.stan',
                             model_str == 'density_dependent' ~ 'mrna_transport_density_dependent_with_blocking.stan'
                            )
    # }
      print(paste('Using the following stan file: ', stan_file, sep=''))
      stan_list = list(y = exp_data,
                       T1  = nSamples,
                       T2 = nSamples+nTest+nTestOE,
                       T3 = nTestOE,
                       y0 = m0,
                       t0 = times$t0$estimate[3],
                       ts1 = times$ts1,
                       ts2 = times$ts2,
                       ts3 = times$ts4,
                       OE_blocked = blocked,
                       OE_producers = producers,
                       y_OE = overexpression_data
                       )
      initF <- function() list(a=10, b=0.2, sigma=rep(1.25,1), nu=0.9, phi=0.35, beta=0.01) #initialising can help speed up the warm up phase   

    estimates <- stan(file = stan_file,
                      data = stan_list,
                      seed = 42,
                      chains = 4,
                      warmup = 1000,
                      iter = 2000,
                      init = initF
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
  
  #look at posterior predictive distn
  source('post_pred_plot.R')
    p1 <- post_pred_plot(test_data,times$ts2,nTest+nSamples+nTestOE,'y_pred',
                         estimates,identifier,title_stem='plots/posterior_pred',
                         ts_test=times$ts3,OE_test=times$ts4,filter_out_wt=TRUE)
    p2 <- post_pred_plot(overexpression_data,times$ts4,nTestOE,'y_pred_OE',
                         estimates,identifier,title_stem='plots/posterior_pred_OE',
                         ts_test=times$ts3,OE_test=times$ts4,filter_out_wt=FALSE)
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
    mcmc_pairs(draws,pars=c('a','b','nu','sigma','phi'), off_diag_fun="hex")
    ggsave(paste('plots/pairs_bayes_',identifier, '.eps',sep=''),device=cairo_ps)
    color_scheme_set("brightblue")
    mcmc_scatter(draws,pars=c('a','b'))
    ggsave(paste('plots/scatter',identifier, '.eps',sep=''),device=cairo_ps)
    
    cairo_ps(paste('plots/pairs',identifier, '.eps',sep=''))
    pairs(estimates, pars = c('a','b','nu','phi','sigma'))
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
