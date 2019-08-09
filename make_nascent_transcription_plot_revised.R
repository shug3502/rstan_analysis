nSamples = 16
nTestOE = 9
library(dplyr)
source('extract_times_and_scaling.R')
times = extract_times_and_scaling(nSamples,1,nTestOE)
source('get_nascent_transcription.R')
##############
#for OE
nascent_transcription = get_nascent_transcription(nTestOE,'OE')
nt_df_oe <- data.frame(nascent_transcription) 
colnames(nt_df_oe) = seq_len(16)
nt_df_oe <- nt_df_oe %>%
  mutate(time = times$ts4) %>%
  gather('cellID','transcription',-time) %>%
  mutate(cellID=as.integer(cellID))

measurements <- read.csv('~/Documents/DPhil/OXMATHSDOCS/FISH_data/rstan_analysis/data/Overexpression5/area_of_transcription_sites.csv')
num_pixels = sum(measurements$Area) / 0.1395^2 #divide by area of x-y pixel

nt_df_oe <- nt_df_oe %>%
  mutate(intensity_per_pixel = transcription/num_pixels)

library(ggplot2)
ggplot(nt_df_oe,aes(x=cellID,y=intensity_per_pixel)) + 
  geom_line() + 
  geom_point() + 
  labs(y = expression(paste("Intensity of transcription foci")), x = 'CellID') +
#  scale_colour_continuous(name="Time\n(hrs)") +
  facet_wrap(~round(time,digits=1)) +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 16),
        strip.text = element_text(size = 8))
ggsave('plots/nascent_trancription_plot_oe_revised.eps')      

####################
#for WT
nascent_transcription = get_nascent_transcription(nSamples,'WT')
area_in_nucleus = 1.2
num_pixels = area_in_nucleus / 0.1395^2 #divide by area of x-y pixel
nt_df_wt <- data.frame(nascent_transcription[7:11,]) 
colnames(nt_df_wt) = seq_len(16)
nt_df_wt <- nt_df_wt %>%
  mutate(time = times$ts1[7:11]) %>%
  gather('cellID','transcription',-time) %>%
  mutate(cellID=as.integer(cellID)) %>%
  mutate(intensity_per_pixel = transcription/num_pixels)

ggplot(nt_df_wt,aes(x=cellID,y=intensity_per_pixel)) + 
  geom_line() + 
  geom_point() + 
  labs(y = expression(paste("Intensity of transcription foci")), x = 'CellID') +
  #  scale_colour_continuous(name="Time\n(hrs)") +
  facet_wrap(~round(time,digits=1)) +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 16),
        strip.text = element_text(size = 8))
ggsave('plots/nascent_trancription_plot_wt.eps')      