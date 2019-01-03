#compare all models with homogeneous overdispersion and inhomogeneous overdispersion
#JH 19/12/18
#######
source('forward_simulate.R')
#######
# inhomogeneous_weights_df <- inhomogeneous_weights_df %>%
#   mutate(overdispersion='inhomogeneous')
# homogeneous_weights_df <- homogeneous_weights_df %>%
#   mutate(overdispersion='homogeneous')
# weights_df <- full_join(inhomogeneous_weights_df,homogeneous_weights_df) %>%
#   mutate(overdispersion=factor(overdispersion))
#######
i<-0; j<-0
aM_inhomogeneous <- list()
simple_models <- models_to_sim[!(models_to_sim %in% models_with_density_dependence)]
expose_stan_functions('mrna_transport_with_blocking.stan')
aM_inhomogeneous[simple_models] <- purrr::map(simple_models, function(mID) alternative_analysis(mID,drop_egg_ind = i,drop_cell_ind = j))
expose_stan_functions('mrna_transport_density_dependent_with_blocking.stan') #should replace previous functions and work with previous wrappers
DD_models <- models_to_sim[models_to_sim %in% models_with_density_dependence]
aM_inhomogeneous[DD_models] <- purrr::map(DD_models, function(mID) alternative_analysis(mID,drop_egg_ind = i,drop_cell_ind = j))
############
source('all_the_same.R')
aM_homogeneous <- list()
simple_models <- models_to_sim[!(models_to_sim %in% models_with_density_dependence)]
expose_stan_functions('mrna_transport_with_blocking.stan')
aM_homogeneous[simple_models] <- purrr::map(simple_models, function(mID) alternative_analysis(mID,drop_egg_ind = i,drop_cell_ind = j))
expose_stan_functions('mrna_transport_density_dependent_with_blocking.stan') #should replace previous functions and work with previous wrappers
DD_models <- models_to_sim[models_to_sim %in% models_with_density_dependence]
aM_homogeneous[DD_models] <- purrr::map(DD_models, function(mID) alternative_analysis(mID,drop_egg_ind = i,drop_cell_ind = j))
############

#now combine
aM_all <- c(aM_homogeneous,aM_inhomogeneous)
#get nicely formatted output from comparison
w_stack = loo::loo_model_weights(aM_all,method='stacking')
w_pseudo_bma = loo::loo_model_weights(aM_all,method="pseudobma")
weights_df <- data_frame(model=rep(models_to_sim,2),
                         pseudo_bma=w_pseudo_bma,stacking=w_stack,
                         overdispersion=c(rep('homogeneous',8),rep('inhomogeneous',8)))

fontsize = 18
method_names <- c(
  `pseudo_bma` = "Pseudo BMA+",
  `stacking` = "Stacking",
  `homogeneous` = "Homogeneous",
  `inhomogeneous` = "Inhomogeneous"  
)

weights_df %>%
  gather(value=weight,key=method,-model, -overdispersion) %>%
  ggplot(aes(x=model,y=weight)) +
  geom_col() +
  facet_grid(rows=vars(method),cols=vars(overdispersion),
             labeller = as_labeller(method_names)) + 
  theme_bw() + 
  theme(text = element_text(size = fontsize), axis.text = element_text(size = fontsize),
        legend.position = "none", strip.text = element_text(size = fontsize)) +
  xlab('Model') + 
  ylab('Weight') + 
  scale_y_continuous(limits=c(0,1))

  ggsave(paste('plots/weights_model_comparison_homogeneous_and_inhomogeneous.eps',sep=''),device=cairo_ps)

