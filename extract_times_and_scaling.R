extract_times_and_scaling <- function(nSamples,nTest,nTestOE,optional_plot=FALSE,test_on_mutant_data=FALSE){
  #Get time series for egg chamber times and fit linear model 
  #Last Edited: 27/06/2018
  #JH
  #################################
  require(dplyr)
  if (test_on_mutant_data){ warning('test_on_mutant_data is an argument no longer used by extract_times_and_scaling')}

  egg_chamber_areas <- rep(0,nSamples+nTest+nTestOE)
  stages <- rep(0,nSamples+nTest+nTestOE) #don't need to extract estimated stages for each egg chamber example
  for (j in seq_len(nSamples+nTest)){
    egg_chamber_areas[j] <- as.numeric(read.table(paste('data/Example',j,'/area.txt',sep='')))
    temp = read.table(file = paste('data/Example',j,'/','filenames.txt',sep='')) %>%
      unlist %>%
      stringr::str_extract(., 'stg.') %>%
      stringr::str_split(.,'stg',simplify=TRUE)
    stages[j] = temp[1,2] %>% as.numeric
  }

  for (j in seq_len(nTestOE)){
    egg_chamber_areas[nSamples+nTest+j] <- as.numeric(read.table(paste('data/Overexpression',j,'/area.txt',sep='')))
    temp = read.table(file = paste('data/Overexpression',j,'/','filenames.txt',sep='')) %>%
      unlist %>%
      stringr::str_extract(., 'stg.') %>%
      stringr::str_split(.,'stg',simplify=TRUE)
    stages[nSamples+nTest+j] = temp[1,2] %>% as.numeric
  }
    
  #take l0 = log(20) as initial time (when using length)
  #take measure egg chamber lengths as a scaled time variable
  log_areas = sort.int(log(egg_chamber_areas[1:nSamples]),index.return=TRUE)
  ts1 = log_areas$x #ts needs to be an ordered time vector
  #log_areas_test = sort.int(log(egg_chamber_areas[(nSamples+1):(nSamples+nTest)]),index.return=TRUE)
  #ts2 = log_areas_test$x #includes the extra test sets or ts = setdiff(ts2,ts1)
  log_areas_test = sort.int(log(egg_chamber_areas),index.return=TRUE)
  ts2 = log_areas_test$x #includes the extra test sets or ts = setdiff(ts2,ts1)
  sort_indices1 = log_areas$ix
  sort_indices2 = log_areas_test$ix
  ts3 = setdiff(ts2,ts1) 
  log_areas3 = sort.int(log(egg_chamber_areas[(nSamples+1):(nSamples+nTest)]),index.return=TRUE)
  sort_indices3 = log_areas3$ix
  log_areas4 = sort.int(log(egg_chamber_areas[(nSamples+nTest+1):(nSamples+nTest+nTestOE)]),index.return=TRUE)
  sort_indices4 = log_areas4$ix
  ts4 = log_areas4$x

  #######################################
  #fit linear models
  lm_time <- lm(log(egg_chamber_areas[1:nSamples]) ~ stages[1:nSamples])
  df <- data.frame(la = log(egg_chamber_areas[1:nSamples]), stages = stages[1:nSamples])
  df <- df %>% mutate(pred_age = predict(lm_time,df))
  real_time <- c(8,6,5,3,6,6) #time for each stage from start of 3 to end of 8, according to jia 2016
  cum_time <- cumsum(real_time)
  ct2 <- (cum_time[1:5]+cum_time[2:6])/2
  df <- df %>% mutate(time_hrs = ct2[stages-3])
  #fit another model to time in hours
  lm_time_hrs = lm(la ~ time_hrs,data=df)
  df <- df %>% mutate(pred_age_hrs = predict(lm_time_hrs,df))      
##########################################  
  if (optional_plot){
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
    ggsave('plots/timescale_model.eps')
  }
 #########################################
##use coefficients of linear model
#   t0 = coef(lm_time)[1] #log(400)
#   time_scaling = coef(lm_time)[2]
##or 
  t0 = coef(lm_time_hrs)[1] #log(400)
  time_scaling = coef(lm_time_hrs)[2]
  return(list(t0=t0,ts1=ts1,ts2=ts2,ts3=ts3,ts4=ts4,sort_indices1=sort_indices1,sort_indices2=sort_indices2,sort_indices3=sort_indices3,sort_indices4=sort_indices4,time_scaling=time_scaling))
}
