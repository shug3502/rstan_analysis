#run model comparison
#JH updated 06/11/18
#######
run_mcmc=TRUE
id = 'v510_minimal'

source('mrna_transport_full.R')
#simple model
identifier <- paste(id,'M0_simple',sep='')
res_M0 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, train_on_OE = FALSE,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
				       model_str='simple')

# decay model
identifier <- paste(id,'decay',sep='')
res_decay = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, train_on_OE = FALSE,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       model_str='decay')

#density dependent model, simple version
identifier <- paste(id,'M4_density_dependent',sep='')
res_M4 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = run_mcmc,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = TRUE,
                                       show_diagnostic_plots = FALSE, train_on_OE = FALSE,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       model_str='density_dependent')

