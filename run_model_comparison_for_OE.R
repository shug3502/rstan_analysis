#run model comparison
#######
run_mcmc=TRUE
omit_OE_data_pts = -1
source('mrna_transport_full.R')
#M0
identifier <- 'v451WT_simple'
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
identifier <- 'v452_with_blocking'
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
identifier <- 'v453_with_production'
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
identifier <- 'v454_with_production4a'
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
identifier <- 'v455WT_with_production4a_a'
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
identifier <- 'v457WT_blocking_and_production'
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
identifier <- 'v458WT_density_dependent_2a'
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
identifier <- 'v459WT_density_dependent_4a'
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
identifier <- 'v460WT_density_dependent_with_blocking'
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
