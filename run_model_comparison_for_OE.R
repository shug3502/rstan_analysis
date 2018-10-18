#run model comparison
#######
run_mcmc=FALSE
omit_OE_data_pts = FALSE #c(-1,-2)
id = 'v470_with_outliers_'  #'v480_without_outliers2_'

source('mrna_transport_full.R')
#M0
identifier <- paste(id,'M0_simple',sep='')
res_M0 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE,
                                       use_binary_producers = FALSE, omit_OE_data_pts=omit_OE_data_pts,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
#M1
identifier <- paste(id,'M1_with_blocking',sep='')
res_M1 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE, omit_OE_data_pts=omit_OE_data_pts,
                                       use_binary_producers = FALSE, use_blocked_RCs = TRUE,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
#M2
identifier <- paste(id,'M2_with_production',sep='')
res_M2 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE,
                                       use_binary_producers = TRUE, omit_OE_data_pts=omit_OE_data_pts,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
#M3
identifier <- paste(id,'M3_with_production2a_a',sep='')
res_M3 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE,
                                       use_binary_producers = 2, omit_OE_data_pts=omit_OE_data_pts,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
#M4
identifier <- paste(id,'M4_with_production_4a_a',sep='')
res_M4 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE,
                                       use_binary_producers = c(4,1), omit_OE_data_pts=omit_OE_data_pts,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
# #M5
# identifier <- 'v446WT_with_production_fit_ka'
# res_M5 = mrna_transport_inference_full(identifier = identifier,
#                                        use_real_data = TRUE, run_mcmc = run_mcmc,
#                                        nSamples = 9, nTest = 11, nTestOE = 9,
#                                        verbose = FALSE, compare_via_loo = TRUE,
#                                        show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
#                                        use_prior_predictive = FALSE, train_on_OE = FALSE,
#                                        use_binary_producers = c(-1,1), omit_OE_data_pts=omit_OE_data_pts,
#                                        parametersToPlot = c('a','b','nu','phi','sigma','alpha'),
#                                        is_nu_uniform = TRUE, no_decay_model = TRUE)
#M6
identifier <- paste(id,'M6_blocking_and_production',sep='')
res_M6 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE, omit_OE_data_pts=omit_OE_data_pts,
                                       use_binary_producers = c(4,1), use_blocked_RCs = TRUE,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
#M7
identifier <- paste(id,'M7_density_dependent_2a',sep='')
res_M7 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE, omit_OE_data_pts=omit_OE_data_pts,
                                       use_binary_producers = c(2,2), use_density_dependence=TRUE,
                                       parametersToPlot = c('a','b','nu','phi','sigma','beta'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
#M8
identifier <- paste(id,'M8_density_dependent_4a',sep='')
res_M8 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE, omit_OE_data_pts=omit_OE_data_pts,
                                       use_binary_producers = c(4,1), use_density_dependence=TRUE,
                                       parametersToPlot = c('a','b','nu','phi','sigma','beta'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
#M9
identifier <- paste(id,'M9_density_dependent_with_blocking',sep='')
res_M9 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = TRUE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE,
                                       use_binary_producers = FALSE, use_blocked_RCs = TRUE,
                                       use_density_dependence=TRUE, omit_OE_data_pts=omit_OE_data_pts,
                                       parametersToPlot = c('a','b','nu','phi','sigma','beta'),
                                       is_nu_uniform = TRUE, no_decay_model = TRUE)
