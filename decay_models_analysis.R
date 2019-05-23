library(loo)
library(rstan)
options(mc.cores = 4)
identifier <- 'v601WTdecay'
source('mrna_transport_full.R')
res_WT = mrna_transport_inference_full(identifier,
                                       use_real_data = TRUE, run_mcmc = FALSE,
                                       nSamples = 9, nTest = 11, nTestOE = 14,
                                       verbose = TRUE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = TRUE, train_on_OE = FALSE,
                                       parametersToPlot = c('a','b','nu','phi','sigma','gamma'),
                                       alpha=2.18, model_str='decay')

loo_1 <- res_WT[[2]]
##################
#and the same for the OE
identifier <- 'v602OEdecay'
res_OE = mrna_transport_inference_full(identifier,
                                       use_real_data = TRUE, run_mcmc = FALSE,
                                       nSamples = 9, nTest = 11, nTestOE = 14,
                                       verbose = TRUE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = TRUE, train_on_OE = TRUE,
                                       parametersToPlot = c('a','b','nu','phi','sigma','gamma'),
                                       alpha=2.18, model_str='decay')

loo_2 <- res_OE[[2]]
#cant compare like this as not using the same data 
#compare(loo_1,loo_2)
