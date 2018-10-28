#see extract_times_and_scaling
get_just_areas <- function(nSamples,nTest,nTestOE){
 egg_chamber_areas <- rep(0,nSamples+nTest+nTestOE)
  for (j in seq_len(nSamples+nTest)){
    egg_chamber_areas[j] <- as.numeric(read.table(paste('data/Example',j,'/area.txt',sep='')))
  }

  for (j in seq_len(nTestOE)){
    egg_chamber_areas[nSamples+nTest+j] <- as.numeric(read.table(paste('data/Overexpression',j,'/area.txt',sep='')))
  }
  return(egg_chamber_areas)
}


get_just_stages <- function(nSamples,nTest,nTestOE){
  stages <- rep(0,nSamples+nTest+nTestOE) #don't need to extract estimated stages for each egg chamber example
  for (j in seq_len(nSamples+nTest)){
    temp = read.table(file = paste('data/Example',j,'/','filenames.txt',sep='')) %>%
      unlist %>%
      stringr::str_extract(., 'stg.') %>%
      stringr::str_split(.,'stg',simplify=TRUE)
    stages[j] = temp[1,2] %>% as.numeric
  }

  for (j in seq_len(nTestOE)){
    temp = read.table(file = paste('data/Overexpression',j,'/','filenames.txt',sep='')) %>%
      unlist %>%
      stringr::str_extract(., 'stg.') %>%
      stringr::str_split(.,'stg',simplify=TRUE)
    stages[nSamples+nTest+j] = temp[1,2] %>% as.numeric
  }
return(stages)
}

convert_to_hrs <- function(stages){
  real_time <- c(8,6,5,3,6,6) #time for each stage from start of 3 to end of 8, according to jia 2016
  cum_time <- cumsum(real_time)
  ct2 <- (cum_time[1:5]+cum_time[2:6])/2 #take midpoint of duration of each stage, as a single stage can be long and do not know where each example falls within this
  cum_time_hrs <- ct2[stages-3]
  return(cum_time_hrs)
}
