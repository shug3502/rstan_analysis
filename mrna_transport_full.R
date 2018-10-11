#setwd('~/Documents/FISH_data/rstan_analysis')
mrna_transport_inference_full <- function(identifier='full_v099',use_real_data=FALSE,run_mcmc=FALSE,nSamples=15,nTest=5,nTestOE=3,
                                             parametersToPlot = c("theta","phi","sigma","a","b"),verbose=FALSE,compare_via_loo=FALSE,
                                             show_diagnostic_plots=FALSE, use_hierarchical_model=FALSE, use_prior_predictive=TRUE,
                                             use_binary_producers=FALSE, use_blocked_RCs=FALSE, train_on_OE=FALSE,
                                             is_nu_uniform=TRUE, no_decay_model=FALSE){
  library(rstan)
  library(mvtnorm)
  library(dplyr)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  if (use_hierarchical_model & (length(which(parametersToPlot=='a'))==0 | length(which(parametersToPlot=='b'))==0)){
    parametersToPlot[which(parametersToPlot=='a')] = 'mu' #replace these for hierarchical model
    parametersToPlot[which(parametersToPlot=='b')] = 'mu'
  } 
  
  #############################################################
  #Fit linear model to log(length) of egg chambers to get time of development.
  #Could use area or volume of egg chamber. Compare to Jia 2016 Sci Reports.
  #Area may be better as in some examples egg chambers are quite squashed.
  #To get lengths read length.txt file in each folder
  
  source('extract_times_and_scaling.R')
  times = extract_times_and_scaling(nSamples,nTest,nTestOE)
  m0 = c(0, rep(0,15)) #initial condition
  #############################################################
  
  #############################################################
  # source('get_nc_transition_matrix.R')
  # B1 = get_nc_transition_matrix(0) %>% as.vector
  # B2 = get_nc_transition_matrix(1) %>% as.vector
  
  #############################################################
  
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
      if (verbose){
        print('swapped')
        print(times$ts1)
        print(nSamples)
      }
    }
  } else {
    warning('TODO: update simulated data for full model')
    #sample from the model to get fake data 
    mc <- stan_model('model_comparison_at_stst_with_decay.stan')  #('mrna_transport5.stan')
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
    #plot(s[1,1,seq(from=1, to=(nSamples*16-1), by=nSamples)])
    #plot(s[1,1,seq(from=nSamples, to=(nSamples*16), by=nSamples)])
    boxplot(s[,1,seq(from=nSamples, to=((nSamples+nTest+nTestOE)*16), by=nSamples)])
    test_data = matrix(s[1,1,1:(16*(nSamples+nTest+nTestOE))],nrow=(nSamples+nTest+nTestOE),byrow=FALSE) #this is our fake data
    exp_data = test_data[!(times$ts2 %in% times$ts3),] 
    overexpression_data = test_data[(times$ts2 %in% times$ts4),]
  }
  if (use_binary_producers){
    source('get_producers.R') #use heterogeneous production information
    producers = 2*get_producers(nTestOE)[times$sort_indices4,] #provides matrix of heterogeneous production due to patch overexpression mutant
  } else {
    producers = matrix(rep(2,nTestOE*16),ncol=16)
    producers[,1] = 0
  }
  print(producers)  
  if (use_blocked_RCs){
    source('get_blocked_indices.R')
    blocked = get_blocked_indices(nTestOE)[times$sort_indices4,]
  } else {
    blocked = matrix(rep(1,nTestOE*3),ncol=3)
  }
  print(blocked)
  
  ##########################
  if (run_mcmc) {
    if (!is_nu_uniform & use_hierarchical_model) warning('Not yet implemented a non uniform hierarchical model. Using normal non-uniform model')
    if (use_prior_predictive){
      stan_file = case_when(
        use_blocked_RCs ~ 'prior_predictive_with_blocking.stan',
        no_decay_model ~ 'prior_predictive_no_decay.stan',
        !is_nu_uniform ~ 'prior_predictive_nu_varying_spatially.stan',
        use_hierarchical_model ~ 'prior_predictive_hierarchical.stan',
        TRUE ~ 'prior_predictive_full.stan')
    } else {
      stan_file = case_when( 
        use_blocked_RCs ~ 'mrna_transport_with_blocking.stan',
        no_decay_model ~ 'mrna_transport_no_decay.stan',
        !is_nu_uniform ~ 'mrna_transport_full_nu_varying_spatially.stan',
        use_hierarchical_model ~ 'mrna_transport_full_hierarchical.stan',
        TRUE ~ 'mrna_transport_full.stan')   #'mrna_transport_reparametrised.stan')
      }
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
    if (!use_hierarchical_model){
      #draw initial values from the prior
      initF <- function() {
	b=abs(10*rnorm(1))
	return(list(a=abs(10*rnorm(1)),b=b,sigma=abs(10*rnorm(1)),nu=runif(1),phi=0.345,gamma=abs(rnorm(1)*b/10)))
      }
      #initF <- function() list(a=9, b=0.18, sigma=1.25, nu=0.9, phi=0.3)    
    } else {
      initF <- function() list(mu=c(9, 0.18, 2.2), sigma=1.25, phi=0.3)    
    }
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
  #pairs(estimates, pars = parametersToPlot)
  #dev.off()
  
  #look at posterior predictive distn
  source('post_pred_plot.R')
  source('post_pred_animation.R')
  if (use_prior_predictive) {
    p1 <- post_pred_plot(exp_data, times$ts1, nSamples, 'y_pred', estimates, identifier, title_stem='plots/prior_pred')
  } else {
    p1 <- post_pred_plot(test_data,times$ts2,nTest+nSamples+nTestOE,'y_pred',
                         estimates,identifier,title_stem='plots/posterior_pred',
                         ts_test=times$ts3,OE_test=times$ts4,filter_out_wt=TRUE)
    p2 <- post_pred_plot(overexpression_data,times$ts4,nTestOE,'y_pred_OE',
                         estimates,identifier,title_stem='plots/posterior_pred_OE',
                         ts_test=times$ts3,OE_test=times$ts4,filter_out_wt=FALSE)
    # p3 <- post_pred_animation(overexpression_data,times$ts4,nTestOE,'y_pred_OE',
    #                      estimates,identifier,title_stem='plots/posterior_pred_OE',
    #                      ts_test=times$ts3,OE_test=times$ts4)
    # p4 <- post_pred_animation(test_data,times$ts2,nTest+nSamples+nTestOE,'y_pred',
    #                      estimates,identifier,title_stem='plots/posterior_pred',
    #                      ts_test=times$ts3,OE_test=times$ts4)
  }  
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
