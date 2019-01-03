source('get_just_areas.R')
wt_fluorescence <- read.csv('data/total_fluorescence_wt.csv')
oe_fluorescence <- read.csv('data/total_fluorescence_oe.csv')

wt_total <- wt_fluorescence$RawIntDen
oe_total <- oe_fluorescence$RawIntDen
nSamples <- nrow(wt_fluorescence) 
nTestOE <- nrow(oe_fluorescence)
wt_log_areas <- get_just_areas(nSamples,0,0)
oe_log_areas <- get_just_areas(0,0,nTestOE)

 lm(wt_total ~ wt_log_areas)
 lm(oe_total ~ oe_log_areas)

library(ggplot2)
library(dplyr)
df <- data_frame(phenotype = c(rep('wt',nSamples),rep('oe',nTestOE)),
                 fluorescence = c(wt_total,oe_total),
                 log_areas = c(wt_log_areas,oe_log_areas))
ggplot(df,aes(y=fluorescence,x=log_areas,group=phenotype,color=phenotype)) + 
  geom_smooth(method=lm) + 
  geom_point()

#The gradients of the lines seem to be different, with the gradient of the overexpression line approximately half of that of the wild type.
#Are there other imaging settings such as power that are not the same between phenotypes?
#
