post_pred_plot <- function(raw_data,ts,nSamples,params,estimates,identifier,
                           title_stem='plots/posterior_pred',ts_test=vector(),
                           OE_test=vector(),neighbouring_ncs_only=FALSE,
                           filter_out_wt=FALSE){
require(tidyr)
require(ggplot2)
require(magrittr)
pred <- as.data.frame(estimates, pars = params) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.025),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.975))
xdata <- data.frame(rna = as.vector(raw_data),cellID = as.vector(matrix(rep(1:16,nSamples),nrow=nSamples,byrow=TRUE)),time = rep(ts,16))
pred %<>% bind_cols(xdata) %>%
        mutate(split = case_when(time %in% OE_test ~ 'overexpression',
                                 time %in% ts_test ~ 'test',
                                 TRUE ~ 'train')) 
if (filter_out_wt){
pred %<>% filter(split != 'overexpression') %>%
        filter(!(split == 'test' & rna==0))
}

if (neighbouring_ncs_only) pred <- pred %>% filter(cellID %in% c(1,2,3,5,9)) #only want to plot nieghbouring cells to oocyte
p1 <- ggplot(pred, aes(x = time, y = median)) +
  geom_line() +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
  facet_wrap(~cellID,scales='free') +   #needed to remove factor(cellID) 
  labs(x = "Time (hrs)", y = "mRNA") +
  theme_bw() +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))
if (!is.na(raw_data[1])){
p1 <- p1 + geom_point(aes(x = time, y = rna, colour = split))
}
print(p1)
ggsave(paste(title_stem,identifier, '.eps',sep=''),device=cairo_ps)
return(p1)
}
