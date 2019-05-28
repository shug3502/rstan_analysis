# make key graph
# JH updated 06/11/18
library(dplyr)
library(ggplot2)
library(tidyr)
nTestOE=14
source('extract_times_and_scaling.R')
times = extract_times_and_scaling(16,4,nTestOE)
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
font_size = 32
h1 <- ggplot(data = quantify_for_ilan, aes(x=time,y=number,color=phenotype)) + 
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  theme(text = element_text(size = font_size), axis.text = element_text(size = font_size),
        strip.text = element_text(size = font_size)) +
  xlab('Time') + 
  ylab('mRNA complexes\nin oocyte')
print(h1)
ggsave(paste('plots/number_in_oocyte_','.eps',sep=''),device=cairo_ps, height=9,width=12)
