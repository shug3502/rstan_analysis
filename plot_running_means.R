
plot_running_means <- function(estimates,parametersToPlot=c("phi","nu"),identifier=""){
  #make additional running mean plots from mcmc output
  #estimates should be a stanfit object
  library(ggmcmc)
  library(rstan)
#  estimates <- readRDS(paste('fits/mrna_transport_estimates',identifier,'.rds',sep=''))
  S <- ggs(rstan::As.mcmc.list(estimates,pars=parametersToPlot))
  p1 <- ggs_running(S,greek=TRUE) +
    theme_bw() + 
    scale_x_continuous(breaks = c(0,1000)) +
    theme(axis.text.x = element_text(angle = 90, hjust = -1))
  print(p1)
ggsave(paste("plots/running_mean_",identifier,".eps",sep=""))

p2 <- ggs_traceplot(S,greek=TRUE) + 
  theme_bw() + 
  scale_x_continuous(breaks=c(0,1000))
print(p2)
ggsave(paste("plots/traceplot_",identifier,".eps",sep=""))
}
