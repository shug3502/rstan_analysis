#use the quantitative observed data on nascent transcription to estimate transcription rates for the overexpressor data
#note our observations only capture a certain point in time
#we assume all the data is drawn from some characteristic distribution of transcription levels across the cells
############################

estimate_adjusted_producers <- function(nTestOE) {
  library(rstan)
  library(dplyr)
  library(tidyr)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  source('extract_times_and_scaling.R')
  times = extract_times_and_scaling(1,1,nTestOE)
  source('get_nascent_transcription.R')
  #for OE
  nascent_transcription = get_nascent_transcription(nTestOE,'OE')
  nt_df_oe <- data.frame(nascent_transcription) 
  colnames(nt_df_oe) = seq_len(16)
  nt_df_oe <- nt_df_oe %>%
    mutate(time = times$ts4) %>%
    gather('cellID','transcription',-time) %>%
    mutate(cellID=as.integer(cellID))
    
  obs_transcription <- nt_df_oe %>%
    select(cellID,transcription,time) %>%
    spread(cellID,transcription) %>%
    ungroup() %>%
    select(-time) %>% 
    as.matrix()
  stan_file <- "missing_data_poisson_transcription.stan"
  stan_list <- list(T=nTestOE,
                    transcription_obs=obs_transcription[,2:16]
  )
  estimates <- stan(file = stan_file,
                    data = stan_list,
                    seed = 42,
                    chains = 4,
                    warmup = 1000,
                    iter = 2000
  )
  estimated_producers = rstan::extract(estimates,pars='a',permuted=TRUE)[[1]] %>% apply(.,2,median) 
  estimated_producers = 2*15*estimated_producers/sum(estimated_producers)
  producers = cbind(rep(0,nTestOE),matrix(rep(estimated_producers,nTestOE),ncol=16))
  return(producers)
}
