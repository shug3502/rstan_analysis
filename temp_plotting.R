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

########
marginal_nu_comparison <- function(){
library(bayesplot)
library(ggridges)
color_scheme_set('purple')
color_scheme_get('purple')
identifier = 'no_decay_v331'
estimates = readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep='')) 
wt.draws <- as.array(estimates, pars='nu')
identifier = 'full_nu_uniform_OE_v302'
estimates = readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep='')) 
oe.draws <- as.array(estimates, pars='nu')
nu_df <- data_frame(wt=as.numeric(wt.draws)[1:2000],oe=as.numeric(oe.draws))
nu_df %>% tidyr::gather(key='phenotype',value='nu') %>%
  ggplot(aes(y=phenotype,x=nu)) +
  geom_density_ridges2(fill='#e5cce5') + 
  theme_bw() +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                    legend.position = "none", strip.text = element_text(size = 12),
                    plot.title = element_text(size=12), 
                    strip.text.x = element_text(size = 16, face="italic")) +
  labs(y='Phenotype',x=expression(paste('Transport bias  ', nu)))
ggsave('plots/ridge_nu_plot.eps',device=cairo_ps)
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
  labs(title='b)')
marginal_WT = marginal_ab_plot()
marginal_WT + res_WT + plot_layout(ncol=1,height=c(1,2))
ggsave(paste('plots/fig5',identifier,'.eps',sep=''),device=cairo_ps)


