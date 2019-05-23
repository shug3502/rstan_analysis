library(dplyr)
library(tidyr)
library(broom)
source('get_just_areas.R')

#manually construct data frame from data in Table S2 in shimada et al 2011
shimada_df = data_frame(stages = 4:8,
                        num_egg_chambers = c(33,42,37,30,20),
                        egg_chamber_areas = c(1912.3, 3399.3, 4576.1, 6096.9, 9884.0), #mean value
                        egg_chamber_area_error = c(65.9, 124.2, 110.3, 291.5, 394.4) #currently not accounting for this
                        )


lm_time <- lm(log(egg_chamber_areas) ~ stages,data=shimada_df)

#df <- shimada_df %>% mutate(pred_age = predict(lm_time,df))
shimada_df <- shimada_df %>% mutate(time_hrs = convert_to_hrs(stages)) #convert to hrs moved to 'get_just_areas.R'
#fit model to time in hours
lm_time_hrs = lm(time_hrs ~ log(egg_chamber_areas),data=shimada_df)
shimada_df <- shimada_df %>% mutate(pred_age_hrs = predict(lm_time_hrs,shimada_df)) 

library(ggplot2)
font_size = 24
g_shimada <- ggplot(data = shimada_df, aes(x=time_hrs,y=log(egg_chamber_areas))) +
  geom_point(shape=2,color='red') + 
  geom_line(color='black',aes(y=log(egg_chamber_areas),x=pred_age_hrs)) +
  theme_bw() + 
  labs(x='age (hrs)',y='log(area)') +
  theme(text = element_text(size = font_size), axis.text = element_text(size = font_size),
        strip.text = element_text(size = 8))


