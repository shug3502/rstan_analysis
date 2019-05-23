
plot_running_means <- function(estimates,parametersToPlot=c("phi","nu")){
  #make additional running mean plots from mcmc output
  #estimates should be a stanfit object
  library(ggmcmc)
  library(rstan)
#  estimates <- readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep=''))
  S <- ggs(rstan::As.mcmc.list(estimates,pars=parametersToPlot))
  p1 <- ggs_running(S,greek=TRUE)
  print(p1)
}
