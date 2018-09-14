library(rstan)
library(tidyr)
mc <- stan_model('model_comparison5.stan')
expose_stan_functions(mc) #get hold of functions defined in the stan code
#keep all other parameters constant during this experiment
m0 = c(0, rep(1,15)) #initial condition
th = c(0.2,10) #c(96.5,249)
sig = 1
phi = 0.3
nSamples = 15
nTest = 5
nTestOE = 0
nTotal = nSamples + nTest
source('extract_times_and_scaling.R')
times = extract_times_and_scaling(nSamples,nTest,nTestOE,test_on_mutant_data=FALSE)
data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
raw_data = data[times$sort_indices2,]
#take a vector of transport bias values for nu to loop through
nu_vec = seq(from=0.7, to=1, by=0.05)
all_extracted_samples = data.frame(rna=numeric(),time=numeric())
for (j in seq_along(nu_vec)){
B = construct_matrix(nu_vec[j])
samples <- stan(file = 'mrna_transport6.stan',
                  data = list (
                    T  = nTotal,
                    y0 = m0,
                    t0 = times$t0$estimate[2],
                    ts = times$ts2,
                    theta = array(th, dim = 2),
                    sigma = sig,
                    phi = phi,
                    B = B %>% as.vector()
                  ),
                  algorithm="Fixed_param",
                  seed = 42,
                  chains = 1,
                  iter =100, 
                  refresh = -1
  )
estimates = rstan::extract(samples, pars='y_ode')
#make the samples go in a 'nice' format
pred <- as.data.frame(estimates, pars = 'y_ode') %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95))
xdata <- data.frame(rna = as.vector(raw_data),cellID = as.vector(matrix(rep(1:16,nTotal),nrow=nTotal,byrow=TRUE)),time = rep(times$ts2,16))
pred <- pred %>% bind_cols(xdata) %>%
        mutate(split = if_else(time %in% times$ts3,'train','test'),
               nu = nu_vec[j])
all_extracted_samples = full_join(all_extracted_samples, pred)
}

p1 <- ggplot(all_extracted_samples, aes(x = time, y = median, group = factor(nu), color=factor(nu)))
p1 <- p1 + geom_line() +
  facet_wrap(~cellID,scales='free') +   #needed to remove factor(cellID) 
  labs(x = "Time (hrs)", y = "mRNA") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        strip.text = element_text(size = 8), legend.position = 'None') + 
  scale_color_discrete(name = "nu") + 
  scale_colour_brewer(palette=7) + 
#  geom_point(aes(x = time, y = rna)) +
  NULL
p1

p2 <- ggplot(all_extracted_samples %>% filter(nu>0.85 & nu<0.95), aes(x = time, y = median, group=factor(nu))) +
  geom_line() +
  facet_wrap(~cellID,scales='free') +   #needed to remove factor(cellID) 
  labs(x = "Time (hrs)", y = "mRNA") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        strip.text = element_text(size = 8), legend.position = 'None') + 
  scale_color_discrete(name = "nu") + 
  scale_colour_brewer(palette=7) + 
#  geom_point(aes(x = time, y = rna)) +
  NULL
p2
ggsave('plots/ODE_sims.eps',device=cairo_ps)
