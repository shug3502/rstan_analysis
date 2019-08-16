#plot unnormalised data distribution

library(dplyr)
library(ggplot2)
library(tidyr)
nTestOE = 14
font_size = 10
source('extract_times_and_scaling.R')
times = extract_times_and_scaling(9,11,nTestOE)
data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
exp_data = data[times$sort_indices1,] #need to sort time series and correspondingly reorder rows
exp_data[is.na(exp_data)]=0 #stan can't deal with NAs
overexpression_data = matrix(as.numeric(read.csv('data/exp_data_overexpression.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
overexpression_data = overexpression_data[times$sort_indices4,] #need to sort time series and correspondingly reorder rows
colnames(exp_data) = seq_len(16)
df_wt <- as.data.frame(exp_data) %>% mutate(time=round(times$ts1,digits=2)) %>% gather(key=cell,value=rna,-time) %>% mutate(phenotype='WT',cell=as.integer(cell))


colnames(overexpression_data) = seq_len(16)
df_oe <- as.data.frame(overexpression_data) %>% mutate(time=round(times$ts4,digits=2)) %>% gather(key=cell,value=rna,-time) %>% mutate(phenotype='OE',cell=as.integer(cell))  

h1 <- full_join(df_wt,df_oe) %>%
  ggplot(aes(x=cell,y=rna,group=time,color=time)) +
  geom_line() +
  facet_grid(~phenotype) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_bw() +
  theme(legend.position='None', strip.text = element_text(size = 10), 
        text = element_text(size = font_size)) + 
  labs(x='Cell ID',y='mRNA')
print(h1)
ggsave('plots/unnormalised_data_distn.eps')

h2 <- full_join(df_wt,df_oe) %>% 
  ggplot(aes(x=cell,y=rna,color=factor(phenotype))) +
  geom_line() +
  facet_wrap(~time) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme_bw() +
  theme(legend.position='None', strip.text = element_text(size = 10), 
        text = element_text(size = font_size)) + 
  labs(x='Cell ID',y='mRNA')
print(h2)
ggsave('plots/unnormalised_data_distn_faceted.eps',width=9,height=9)

