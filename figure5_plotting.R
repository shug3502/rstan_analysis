#########
marginal_ab_plot <- function(identifier = 'no_decay_v331'){
  estimates = readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep='')) 
  library(bayesplot)
  color_scheme_set("purple")
  g1<- mcmc_dens(as.array(estimates),pars=c('a','b')) +
    theme_bw() +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none",
          plot.title = element_text(size=12), 
          strip.text.x = element_text(size = 16, face="italic")) +
    labs(title='a)')
  ggsave(paste('plots/post_marginal_',identifier,'.eps',sep=''),device=cairo_ps)
  return(g1)
}

###################

library(patchwork)
source('mrna_transport_full.R')
identifier = 'no_decay_v331'
res_WT = mrna_transport_inference_full(identifier,
                                       use_real_data = TRUE, run_mcmc = FALSE,
                                       nSamples = 9, nTest = 11, nTestOE = 9,
                                       verbose = FALSE, compare_via_loo = FALSE,
                                       show_diagnostic_plots = FALSE, use_hierarchical_model = FALSE,
                                       use_prior_predictive = FALSE, train_on_OE = FALSE,
                                       is_nu_uniform = TRUE, no_decay_model = TRUE) +
  labs(title='b)') + theme(title=element_text(size=12))
marginal_WT = marginal_ab_plot(identifier)
marginal_WT + res_WT + plot_layout(ncol=1,height=c(1,3))
ggsave(paste('plots/fig5',identifier,'.eps',sep=''),width = 9, height = 9, device=cairo_ps)


