#attempt to do model comparison
#model index is one of the parameters. This has to be a real. 0.5 is the switch between model1 and model2
setwd('~/Documents/FISH_data/rstan_analysis')
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
identifier = 'MCv007' #run identifier
use_real_data <- FALSE

deltaT = 1
nSamples = 10
m0 = c(0, rep(1,15)) #initial condition
th = c(0.2,0.1)
sig = 2.2
t0 = 0.0
ts = seq(deltaT,nSamples * deltaT,deltaT)

#############################################################
#easier to write functions to get transition matrices in R and feed in as data
get_nc_transition_matrix <- function(isBidirectional){
  if (isBidirectional){
    v <- c(c(2,3,5,9),c(1,4,6,10),c(1,7,11),c(2,8,12),c(1,13),c(2,14),c(3,15),c(4,16),1,2,3,4,5,6,7,8)
  } else {
    v <- c(0,1,1,2,1,2,3,4,1,2,3,4,5,6,7,8)
  }
  B = matrix(rep(0,16^2),nrow=16)
  for (j in 1:16){
    B[v[j],j] = 1;
    B[j,j] = -1;
  }
  return(B)
}
B1 = get_nc_transition_matrix(0)
B2 = get_nc_transition_matrix(1)
#############################################################

if (use_real_data){
  #data taken from spot detection on 3 egg chambers
  print('using real data \n')
  exp_data <- matrix(c(708, 47, 117, 50, 91, 34, 86, 49, 76, 94, 100, 63, 110, 73, 81, 71,
                       1805, 283, 194, 168, 382, 179, 67, 22, 85, 70, 56, 47, 174, 136, 20, 17,
                       1454, 261, 144, 180, 108, 142, 10, 96, 18, 127, 20, 53, 19, 109, 20, 11),nrow=nSamples, byrow=TRUE) 
} else {
  #sample from the model to get fake data
  print('using fake generated data')
  samples <- stan(file = 'mrna_transport2.stan',
                  data = list (
                    T  = nSamples,
                    y0 = m0,
                    t0 = t0,
                    ts = ts,
                    theta = array(th, dim = 2),
                    sigma = sig,
                    refresh = -1
                  ),
                  algorithm="Fixed_param",
                  seed = 42,
                  chains = 1,
                  iter =1
  )
  
  s <- rstan::extract(samples,permuted=FALSE)
  plot(s[1,1,seq(from=1, to=(1+(nSamples-1)*16), by=nSamples)])
  plot(s[1,1,seq(from=nSamples, to=(nSamples*16), by=nSamples)])
  exp_data = matrix(s[1,1,1:(16*nSamples)],nrow=nSamples,byrow=FALSE) #this is our fake data
}

##########################

estimates <- stan(file = 'model_comparison_exponentiated.stan',
                  data = list (
                    y = cbind(exp_data,rep(0,nSamples)),
                    T  = nSamples,
                    y0 = c(m0,0),
                    t0 = t0,
                    ts = ts,
                    biB = B1,
                    uniB = B2
                  ),
                  seed = 42,
                  chains = 4,
                  iter = 1000,
                  warmup = 500,
                  control = list(adapt_delta = 0.99)
)

parametersToPlot = c("theta","sigma","lp__")

tryCatch({
  estimates@stanmodel@dso <- new("cxxdso") #seems have to do this to be able to save :S
  saveRDS(estimates, file = paste('fits/model_comparison_estimates',identifier,'.rds',sep='')) # to load use readRDS("fit.rds")
}, warning = function(war) {  
  # warning handler picks up where error was generated
  print(paste("WARNING:  ",war))
}, error = function(err) {  
  # error handler picks up where error was generated
  print(paste("NOT SAVED ERROR:  ",err))
}) # END tryCatch
print(estimates, pars = parametersToPlot)

#pairs(estimates, pars = parametersToPlot)

#######################
#visualisation
source('mcmcDensity.R')
mcmcDensity(estimates, parametersToPlot, byChain = TRUE)
ggsave(paste('plots/denisty',identifier, '.eps',sep=''),device=cairo_ps)

library(bayesplot)
draws <- as.array(estimates, pars=parametersToPlot)
mcmc_trace(draws)
ggsave(paste('plots/trace',identifier, '.eps',sep=''),device=cairo_ps)
mcmc_intervals(draws,pars=c('theta[1]','theta[2]','theta[3]','sigma'))
ggsave(paste('plots/intervals',identifier, '.eps',sep=''),device=cairo_ps)
color_scheme_set("purple")
mcmc_areas(draws,pars=c('theta[1]','theta[2]','theta[3]','sigma'))
ggsave(paste('plots/areas',identifier, '.eps',sep=''),device=cairo_ps)
color_scheme_set("brightblue")
mcmc_scatter(draws,pars=c('theta[1]','theta[3]'))
ggsave(paste('plots/scatter',identifier, '.eps',sep=''),device=cairo_ps)

###############################
#look at posterior predictive distn
xdata <- data.frame(rna = as.vector(cbind(exp_data,rep(0,nSamples))),cellID = as.vector(matrix(rep(1:17,nSamples),nrow=nSamples,byrow=TRUE)),time = rep(ts,17))
pred <- as.data.frame(estimates, pars = "y_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(xdata)

p1 <- ggplot(pred, aes(x = time, y = rna))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "rna") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) + 
  facet_wrap(~factor(cellID))
ggsave(paste('plots/posterior_pred',identifier, '.eps',sep=''),device=cairo_ps)

######################
#cos saving wasn't working
e <- rstan::extract(estimates,pars=parametersToPlot,permuted=TRUE)
save(e,file=paste('fits/alt_save',identifier,'.Rsave',sep=''))
