compare_RC_preds <- function(ts,estimates,identifier,title_stem='plots/compare_RC_preds_'){
  require(tidyr)
  require(dplyr)
  require(ggplot2)
  nSamples <- length(ts)
  get_preds <- function(param) {
    pred <- as.data.frame(estimates, pars = param) %>%
      gather(factor_key = TRUE) %>%
      group_by(key) %>%
      summarize(lb = quantile(value, probs = 0.025),
                median = quantile(value, probs = 0.5),
                ub = quantile(value, probs = 0.975))
    xdata <- data.frame(cellID = as.vector(matrix(rep(1:16,nSamples),nrow=nSamples,byrow=TRUE)),time = rep(ts,16))
    pred <- pred %>% bind_cols(xdata)
    return(pred)
  }
  
  pred <- get_preds('y_pred') %>% mutate(altered='No blocking')
  removed_RC_pred <- get_preds('y_pred_altered') %>% mutate(altered='With blocking')
  preds_compared <- full_join(pred,removed_RC_pred) %>%
    group_by(time) %>% 
    mutate(oocyte_rna = median[1]) %>%
    ungroup() %>%
    mutate(normalised_rna = median/oocyte_rna)

  p1 <- ggplot(preds_compared, aes(x = time, y = median)) +
    geom_line(aes(color=altered)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, color=altered), alpha = 0.25) +
    facet_wrap(~cellID,scales='free') +   #needed to remove factor(cellID) 
    labs(x = "time", y = "rna") +
    theme_bw() +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          #legend.position = "none", 
          strip.text = element_text(size = 8))

  print(p1)
  ggsave(paste(title_stem,identifier, '.eps',sep=''),device=cairo_ps)
  
  bp_colors <- bayesplot::color_scheme_get('purple')
  p2 <- ggplot(preds_compared, aes(x=cellID,y=median,group=time,color=time)) +
    geom_line() +
    scale_colour_gradient(low = bp_colors[[1]], high = bp_colors[[2]]) +
    facet_grid(~altered) + 
    labs(x='Cell ID', y='mRNA',color='Time',title) +
    theme_bw() + 
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          #legend.position = "none", 
          strip.text = element_text(size = 12))
  print(p2)
  ggsave(paste('plots/compare_side_by_side',identifier, '.eps',sep=''),device=cairo_ps)

  p3 <- ggplot(preds_compared, aes(x=cellID,y=normalised_rna,group=time,color=time)) +
    geom_line() +
    scale_colour_gradient(low = bp_colors[[1]], high = bp_colors[[2]]) +
    facet_grid(~altered) + 
    labs(x='Cell ID', y='mRNA',color='Time',title) +
    theme_bw() + 
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          #legend.position = "none", 
          strip.text = element_text(size = 12))
  print(p3)
  return(list(p3,p2,p1,preds_compared))
}
