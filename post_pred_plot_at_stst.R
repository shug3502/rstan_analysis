post_pred_plot_at_stst <- function(raw_data,ts,nSamples,params,estimates,identifier,title_stem='plots/posterior_pred_stst',ts_test=vector()){
  require(tidyr)
  require(ggplot2)
  pred <- as.data.frame(estimates, pars = params) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.05),
              median = quantile(value, probs = 0.5),
              ub = quantile(value, probs = 0.95))
  xdata <- data.frame(rna = as.vector(raw_data),cellID = as.vector(matrix(rep(1:16,nSamples),nrow=nSamples,byrow=TRUE)),time = rep(ts,16))
  pred <- pred %>% bind_cols(xdata) %>%
    mutate(split = if_else(time %in% ts_test,'train','test')) 
  p1 <- ggplot(pred, aes(x = cellID, y = median))
  p1 <- p1 + geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
#    facet_wrap(~factor(cellID)) +  
    labs(x = "Cell ID", y = "rna") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8))
  if (!is.na(raw_data[1])){
    p1 <- p1 + geom_point(aes(x = cellID, y = rna, colour = split))
  }
  ggsave(paste(title_stem,identifier, '.eps',sep=''),device=cairo_ps)
return(p1)
}
