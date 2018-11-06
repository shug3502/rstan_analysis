extract_times_and_scaling <- function(nSamples,nTest,nTestOE,optional_plot=FALSE){
  #Get time series for egg chamber times and fit linear model 
  #Last Edited: 29/10/2018
  #JH
  #################################
  require(dplyr)
  require(broom)
  require(purrr)
  source('get_just_areas.R')

rescale_time <- function(log_area, tau, t0, tol=10^-5){
  t_hat = 1/tau*(log_area - t0)
  return(max(t_hat,tol)) #not sure what to do about the times that turn out negative
}

egg_chamber_areas <- get_just_areas(nSamples,nTest,nTestOE)
stages <- get_just_stages(nSamples,nTest,nTestOE)
  
  #take measure egg chamber lengths as a scaled time variable
if (nSamples>0){ #the if statements deal with having 0 for nSamples etc
  log_areas = sort.int(log(egg_chamber_areas[1:nSamples]),index.return=TRUE)
  ts1 = log_areas$x #ts needs to be an ordered time vector
  sort_indices1 = log_areas$ix
  } else {
  log_areas = numeric(0)
  ts1 = numeric(0)
  sort_indices1 = numeric(0)
}
if (nSamples+nTest+nTestOE>0){
  log_areas_test = sort.int(log(egg_chamber_areas),index.return=TRUE)
  ts2 = log_areas_test$x #includes the extra test sets or ts = setdiff(ts2,ts1)
  sort_indices2 = log_areas_test$ix
} else{
  log_areas = numeric(0)  
  ts2 = numeric(0)
  sort_indices2 = numeric(0)
}
if (nTest>0){
  log_areas3 = sort.int(log(egg_chamber_areas[(nSamples+1):(nSamples+nTest)]),index.return=TRUE)
  sort_indices3 = log_areas3$ix
  ts3 = log_areas3$x
} else {
  log_areas3 = numeric(0)
  sort_indices3 = numeric(0)
  ts3 = numeric(0)
}
if (nTestOE>0){
  log_areas4 = sort.int(log(egg_chamber_areas[(nSamples+nTest+1):(nSamples+nTest+nTestOE)]),index.return=TRUE)
  sort_indices4 = log_areas4$ix
  ts4 = log_areas4$x
} else {
  log_areas4 = numeric(0)
  sort_indices4 = numeric(0)
  ts4 = numeric(0)
}

  #######################################
  #fit linear models
  lm_time <- lm(log(egg_chamber_areas) ~ stages)
  df <- data.frame(la = log(egg_chamber_areas), stages = stages, split = c(rep('train',nSamples),rep('test',nTest),rep('Overexpression',nTestOE)))
  df <- df %>% mutate(pred_age = predict(lm_time,df))
  df <- df %>% mutate(time_hrs = convert_to_hrs(stages)) #convert to hrs moved to 'get_just_areas.R'
  #fit another model to time in hours
  lm_time_hrs = lm(la ~ time_hrs,data=df)
  df <- df %>% mutate(pred_age_hrs = predict(lm_time_hrs,df))      
  # fit linear models grouped by split of data. Flies fed at separate times so grew differently
  lm_by_split <- df %>% group_by(split) %>%  do(tidy(lm(la ~ time_hrs, data = .)))
  lm_by_split <- lm_by_split %>% full_join(df %>%
                               do(tidy(lm(la ~ time_hrs, data = .))) %>%
                               mutate(split='all')
                             )
##########################################  
  if (optional_plot){
    require(ggplot2)
    g <- ggplot(data = df, aes(x=stages,y=la)) +
      geom_point() + 
      geom_line(color='red',aes(x=stages,y=pred_age)) +
      theme_bw() + 
      labs(x='stage',y='log(area)')
    print(g) 
    ###################################
    g <- ggplot(data = df, aes(x=time_hrs,y=la)) +
      geom_point() + 
      geom_line(color='red',aes(x=time_hrs,y=pred_age_hrs)) +
      theme_bw() + 
      labs(x='age (hrs)',y='log(area)')
    print(g) 
    ggsave('plots/timescale_model.eps',device=cairo_ps)
    
    g2 <- ggplot(data=df,aes(x=time_hrs,y=la)) + 
      geom_point() +
      geom_smooth(method='lm') +
      theme_bw() +
      labs(x='time (hrs)',y='log(area)')
    print(g2)
    ggsave('plots/timescale_model2.eps',device=cairo_ps)
    
    h <- ggplot(data=df,aes(x=time_hrs,y=la,color=split)) + 
      geom_point() +
      geom_smooth(method='lm') +
      theme_bw() +
      labs(x='time (hrs)',y='log(area)')
    print(h)
    ggsave('plots/timescale_model_split.eps',device=cairo_ps)
  }
 #########################################
##use coefficients of linear model
  # t0 = coef(lm_time_hrs)[1] #log(400)
  # time_scaling = coef(lm_time_hrs)[2]
  t0 = lm_by_split %>%
    ungroup() %>%
    filter(grepl('Intercept',term)) %>%
    dplyr::select(split,estimate) #can get confused with MASS::select if loaded wrongly
  time_scaling = lm_by_split %>%
    ungroup %>%
    filter(grepl('time_hrs',term)) %>%
    dplyr::select(split,estimate)
  
  ts1 = purrr::map_dbl(ts1,function(x) rescale_time(x,time_scaling %>% filter(grepl('train',split)) %>% select(estimate) %>% .$estimate,t0 %>% filter(grepl('train',split)) %>% select(estimate) %>% .$estimate))
  ts2 = purrr::map_dbl(ts2,function(x) rescale_time(x,time_scaling %>% filter(grepl('train',split)) %>% select(estimate) %>% .$estimate,t0 %>% filter(grepl('train',split)) %>% select(estimate) %>% .$estimate))
  ts3 = purrr::map_dbl(ts3,function(x) rescale_time(x,time_scaling %>% filter(grepl('test',split)) %>% select(estimate) %>% .$estimate,t0 %>% filter(grepl('test',split)) %>% select(estimate) %>% .$estimate))
  ts4 = purrr::map_dbl(ts4,function(x) rescale_time(x,time_scaling %>% filter(grepl('Overexpression',split)) %>% select(estimate) %>% .$estimate,t0 %>% filter(grepl('Overexpression',split)) %>% select(estimate) %>% .$estimate))

  #try using separate models for each split of data
  ts_by_split = sort.int(c(ts1,ts3,ts4),index.return=TRUE)
  ts2 = ts_by_split$x
  sort_indices2 = ts_by_split$ix
  
  return(list(t0=t0,
              ts1=ts1,
              ts2=ts2,
              ts3=ts3, 
              ts4=ts4,
              sort_indices1=sort_indices1,
              sort_indices2=sort_indices2,
              sort_indices3=sort_indices3,
              sort_indices4=sort_indices4,
              time_scaling=time_scaling))
}
