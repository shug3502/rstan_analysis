#This estimates the ratio between production for WT and for the OE mutant
#JH 14/12/2018
#################
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
nTestOE=9
source('extract_times_and_scaling.R')
##################
times = extract_times_and_scaling(9,11,nTestOE)
data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
exp_data = data[times$sort_indices1,] #need to sort time series and correspondingly reorder rows
exp_data[is.na(exp_data)]=0 #stan can't deal with NAs
overexpression_data = matrix(as.numeric(read.csv('data/exp_data_overexpression.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
overexpression_data = overexpression_data[times$sort_indices4,] #need to sort time series and correspondingly reorder rows
colnames(exp_data) = seq_len(16)
df_wt <- as.data.frame(exp_data) %>% mutate(time=times$ts1) %>% gather(key=cell,value=rna,-time) %>% mutate(phenotype='WT',cell=as.integer(cell))
colnames(overexpression_data) = seq_len(16)
df_oe <- as.data.frame(overexpression_data) %>% mutate(time=times$ts4) %>% gather(key=cell,value=rna,-time) %>% mutate(phenotype='OE',cell=as.integer(cell))  

gurken_rna_df <- full_join(df_wt,df_oe)
quantify_for_ilan <- gurken_rna_df %>% 
  group_by(time) %>%
  mutate(phenotype = factor(phenotype),
         number = rna[cell==1],
         fraction = number/sum(rna),
         total = sum(rna)) %>%
  filter(cell==1)

ggplot(quantify_for_ilan, aes(y=total,x=time,color=phenotype,group=phenotype)) + 
  geom_smooth(method = lm) + 
  geom_point()

stan_file = 'how_much_more_produced_in_oe.stan'
stan_list = list(T1=9,
                 T2=nTestOE,
                 x=quantify_for_ilan %>% filter(phenotype=='WT') %>% .$total,
                 y=quantify_for_ilan %>% filter(phenotype=='OE') %>% .$total,
                 ts1=quantify_for_ilan %>% filter(phenotype=='WT') %>% .$time,
                 ts2=quantify_for_ilan %>% filter(phenotype=='OE') %>% .$time
)
estimates <- stan(file = stan_file,
                  data = stan_list,
                  seed = 42,
                  chains = 4,
                  warmup = 1000,
                  iter = 2000
)
print(estimates)

#now plot results
library(bayesplot)
parametersToPlot <- c('gamma','alpha')
draws <- as.array(estimates, pars=parametersToPlot)
mcmc_dens(draws,pars='gamma') +
  labs(y = "Density", x = expression(gamma)) +
  theme_bw() +
  theme(text = element_text(size = 32), axis.text = element_text(size = 32),
        strip.text = element_text(size = 32))
ggsave('plots/production_ratio_gamma.eps',device=cairo_ps)

