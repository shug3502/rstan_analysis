get_adjusted_producers <- function(nTestOE,with_plot=FALSE) {
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(ggplot2)
  source('extract_times_and_scaling.R')
  times = extract_times_and_scaling(1,1,nTestOE)
  source('get_nascent_transcription.R')
  #for OE
  nascent_transcription = get_nascent_transcription(nTestOE,'OE')
  nt_df_oe <- data.frame(nascent_transcription) 
  colnames(nt_df_oe) = seq_len(16)
  nt_df_oe %<>%
    mutate(time = times$ts4) %>%
    gather('cellID','transcription',-time) %>%
    mutate(cellID=as.integer(cellID))
  if (with_plot){
    ggplot(data=nt_df_oe,aes(x=cellID,y=transcription,color=time,group=time)) +
      geom_line() + theme_bw()
  }
  #Now process such that the production data can be used for production models
  nt_df_oe %<>%
    group_by(time) %>% 
    mutate(normalised_transcription = transcription/sum(transcription)*15*2)

  if (with_plot) {
    ggplot(data=nt_df_oe,aes(x=cellID,y=normalised_transcription,color=time,group=time)) +
    geom_line() +
    theme_bw()
  }  
#now convert back to a matrix to put into the model
  adjusted_producers <- nt_df_oe %>%
    select(-transcription) %>%
    spread(cellID,normalised_transcription) %>%
    ungroup() %>%
    select(-time) %>% 
    as.matrix()
}