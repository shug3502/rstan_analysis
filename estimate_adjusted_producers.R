#use the quantitative observed data on nascent transcription to estimate transcription rates for the overexpressor data
#note our observations only capture a certain point in time
#we assume all the data is drawn from some characteristic distribution of transcription levels across the cells
############################

estimate_adjusted_producers <- function(nTestOE,gamma=2,plot_option=F) {
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
                    iter = 2000,
                    refresh=-1
  )
  estimated_producers = rstan::extract(estimates,pars='a',permuted=TRUE)[[1]] %>% apply(.,2,median) 
  estimated_producers = gamma*15*estimated_producers/sum(estimated_producers)
  producers = cbind(rep(0,nTestOE),matrix(rep(estimated_producers,nTestOE),ncol=15,byrow = TRUE))
  
  if (plot_option){
    library(ggplot2)
    ggplot(nt_df_oe,aes(x=cellID,y=transcription/10^5,group=time,color=time)) + 
      geom_line() + 
      geom_point() + 
      labs(y = expression(paste("Trancription (",x, 10^5,')')), x = 'CellID') +
      scale_colour_continuous(name="Time\n(hrs)") +
      theme_bw() +
      theme(text = element_text(size = 32), axis.text = element_text(size = 32),
            strip.text = element_text(size = 32))
    ggsave('plots/nascent_trancription_plot.eps',device=cairo_ps)      

        #also plot the inferred distribution of transcription across cells
    inferred_df = data_frame(a = producers[1,], CellID=seq_len(16))
    ggplot(inferred_df,aes(x=CellID,y=a)) + 
      geom_point(shape=8,size=8) + 
      labs(y = 'Production', x = 'CellID') +
      theme_bw() +
      theme(text = element_text(size = 32), axis.text = element_text(size = 32),
            strip.text = element_text(size = 32))
    ggsave('plots/inhomogeneous_production_inferred.eps',device=cairo_ps)

    #   #now also show the uncertainty in the estimates
  # draws = as.data.frame(rstan::extract(estimates,pars='a',permuted=TRUE)[[1]])
  # names(draws) = seq_len(16)[-1]
  # tidy_draws <- draws %>%
  #   gather('cellID','a') %>%
  #   mutate(cellID = as.integer(cellID))
  # ggplot(tidy_draws,aes(x=factor(cellID),y=a/10^5)) + 
  #   geom_violin()
     } 
  
  return(producers)
}
