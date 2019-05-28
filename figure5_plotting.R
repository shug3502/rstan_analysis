#########
marginal_ab_plot <- function(identifier){
  estimates = readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep='')) 
  library(bayesplot)
  fontsize = 24
  color_scheme_set("purple")
  g1<- mcmc_dens(as.array(estimates),pars=c('a','b')) +
    theme_bw() +
    theme(text = element_text(size = fontsize), axis.text = element_text(size = fontsize),
          legend.position = "none",
          plot.title = element_text(size=fontsize), 
          strip.text.x = element_text(size = fontsize, face="italic")) +
    labs(title='a)')
  ggsave(paste('plots/post_marginal_',identifier,'.eps',sep=''),device=cairo_ps)
  return(g1)
}

###################
plot_figure5 <- function(identifier){
library(ggplot2)
library(patchwork)
source('mrna_transport_full.R')
#id = 'v500_minimal'
#identifier <- paste(id,'M0_simple',sep='')
res_M0 = mrna_transport_inference_full(identifier = identifier,
                                       use_real_data = TRUE, run_mcmc = FALSE,
                                       nSamples = 16, nTest = 4, nTestOE = 14,
                                       verbose = FALSE, compare_via_loo = FALSE,
                                       show_diagnostic_plots = TRUE, train_on_OE = FALSE,
                                       parametersToPlot = c('a','b','nu','phi','sigma'),
                                       model_str='simple')
fontsize = 24
res_WT = res_M0 +
  labs(title='b)') + theme(title=element_text(size=fontsize*2/3))
marginal_WT = marginal_ab_plot(identifier)
marginal_WT + res_WT + plot_layout(ncol=1,height=c(1,3))
ggsave(paste('plots/fig5',identifier,'.eps',sep=''),width = 9, height = 9, device=cairo_ps)

##################
#for pairs plot
library(bayesplot)
bayesplot_theme_set(theme_bw())
bayesplot_theme_update(text = element_text(size = 16, family = "serif"))
estimates = readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep='')) 
color_scheme_set("purple")
draws <- as.array(estimates, pars=c('a','b','nu','sigma','phi'))
attr(draws, "dimnames")[['parameters']] <- c('a','b','ν','σ','φ')
h <- mcmc_pairs(draws,pars=c('a','b','ν','σ','φ'),
                off_diag_fun = 'hex')
print(h)
#easiest to save manually for mcmc_pairs plots, its a grid of ggplot objects in a confusing way
#ggsave(paste('plots/pairs_bayes_',identifier,'.eps',sep=''),device=cairo_ps)
}
################


