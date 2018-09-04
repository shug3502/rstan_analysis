post_pred_animation <- function(raw_data,ts,nSamples,params,estimates,identifier,title_stem='plots/posterior_pred',ts_test=vector(),OE_test=vector(),neighbouring_ncs_only=FALSE){
  library(tidyr)
  library(ggplot2)
  library(gganimate)
  key_col = as.data.frame(estimates, pars = params) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.025),
              median = quantile(value, probs = 0.5),
              ub = quantile(value, probs = 0.975)) %>% 
    select(key)
  xdata <- data.frame(rna = as.vector(raw_data),cellID = as.vector(matrix(rep(1:16,nSamples),nrow=nSamples,byrow=TRUE)),time = rep(ts,16),key=key_col)
  pred <- as.data.frame(estimates, pars = params) %>%
    gather(factor_key = TRUE) %>%
    full_join(xdata) %>%
    mutate(split = case_when(time %in% OE_test ~ 'overexpression',
                             time %in% ts_test ~ 'test',
                             TRUE ~ 'train'))
  
  if (neighbouring_ncs_only) pred <- pred %>% filter(cellID %in% c(1,2,3,5,9)) #only want to plot nieghbouring cells to oocyte
  # p1 <- ggplot(pred, aes(x = cellID, y = median)) +
  #   geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
  #   theme(text = element_text(size = 12), axis.text = element_text(size = 12),
  #         legend.position = "none", strip.text = element_text(size = 8))
  # if (!is.na(raw_data[1])){
  #   p1 <- p1 + geom_point(aes(x = cellID, y = rna, colour = split))
  # }
  # p1 <- p1 + labs(title = 'Time: {frame_time}', x = "cell", y = "rna") +
  #   transition_time(time) +
  #   ease_aes('linear')
  # print(p1)
#  ggsave(paste(title_stem,identifier, '.eps',sep=''),device=cairo_ps)
  return(pred)
}
