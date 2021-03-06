
simulate_from_ODE_model <- function(th = c(0.2,10),
                                    sig = 1,
                                    phi = 0.345,
                                    nSamples = 15,
                                    nTest = 5,
                                    nTestOE = 0){
  library(rstan)
  library(tidyr)
  library(ggimage)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  mc <- stan_model('mrna_transport_with_blocking.stan')
  expose_stan_functions(mc) #get hold of functions defined in the stan code

  y0 = c(0, rep(0,15)) #initial condition
  nTotal = nSamples + nTest
  source('extract_times_and_scaling.R')
  times = extract_times_and_scaling(nSamples,nTest,nTestOE)
  data = matrix(as.numeric(read.csv('data/exp_data.csv',sep=',',header=FALSE,stringsAsFactors = FALSE)),ncol=16,byrow=TRUE)
  raw_data = data[times$sort_indices2,]
  #take a vector of transport bias values for nu to loop through
  nu_vec = seq(from=0.7, to=1, by=0.05)
  all_extracted_samples = data.frame(rna=numeric(),time=numeric())
  for (j in seq_along(nu_vec)){
  B = construct_matrix(nu_vec[j],0)
  samples <- stan(file = 'mrna_transport6.stan',
                      data = list (
                        T  = nTotal,
                        y0 = y0,
                        t0 = times$t0$estimate[2],
                        ts = times$ts2,
                        theta = array(th, dim = 2),
                        sigma = sig,
                        phi = phi,
                        B = B %>% as.vector()
                      ),
                      algorithm="Fixed_param",
                      seed = 42,
                      chains = 1,
                      iter =100, 
                      refresh = -1
      )
    estimates = rstan::extract(samples, pars='y_ode')
    #make the samples go in a 'nice' format
    pred <- as.data.frame(estimates, pars = 'y_ode') %>%
      gather(factor_key = TRUE) %>%
      group_by(key) %>%
      summarize(lb = quantile(value, probs = 0.025),
                median = quantile(value, probs = 0.5),
                ub = quantile(value, probs = 0.975))
    xdata <- data.frame(rna = as.vector(raw_data),cellID = as.vector(matrix(rep(1:16,nTotal),nrow=nTotal,byrow=TRUE)),time = rep(times$ts2,16))
    pred <- pred %>% bind_cols(xdata) %>%
            mutate(split = if_else(time %in% times$ts3,'train','test'),
                   nu = nu_vec[j])
  all_extracted_samples = full_join(all_extracted_samples, pred)
  }
  return(all_extracted_samples)
}

library(ggplot2)
animate_on <- FALSE
my_device <- NULL #or cairo_ps, not sure why sometimes this works or not
font_size <- 14
all_extracted_samples <- simulate_from_ODE_model() 
p1 <- ggplot(all_extracted_samples, aes(x = time, y = median, group = factor(nu), color=factor(nu)))
p1 <- p1 + geom_line() +
  facet_wrap(~cellID,scales='free_y') +   #needed to remove factor(cellID)
  labs(x = "Time (hrs)", y = "mRNA") +
  theme_bw() +
  theme(text = element_text(size = font_size), axis.text = element_text(size = font_size),
        strip.text = element_text(size = 8), legend.position = 'None') +
  scale_color_discrete(name = "nu") +
  scale_colour_brewer(palette=7) +
#  geom_point(aes(x = time, y = rna)) +
  NULL
print(p1)

p2 <- ggplot(all_extracted_samples %>% filter(nu>0.85 & nu<0.95), aes(x = time, y = median, group=factor(nu))) +
  geom_line() +
  facet_wrap(~cellID,scales='free_y') +   #needed to remove factor(cellID) 
  labs(title='a)', x = "Time (hrs)", y = "mRNA\ncomplexes") +
  theme_bw() +
  theme(text = element_text(size = font_size), axis.text = element_text(size = 14),
        strip.text = element_text(size = 8), legend.position = 'None') + 
  scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_color_discrete(name = "nu") + 
  scale_colour_brewer(palette=7) + 
#  geom_point(aes(x = time, y = rna)) +
  NULL
print(p2)
ggsave('plots/fig2a.eps',device = my_device)

p3 <- ggplot(all_extracted_samples %>% filter(nu>0.85 & nu<0.95) %>% filter(time==max(time)), aes(x = cellID, y = median, group=factor(nu))) +
  geom_line() +
  geom_point(shape=8) +
  labs(x = "Cell ID", y = "mRNA\ncomplexes") +
  scale_x_continuous(breaks=c(1,4,8,12,16)) +
  theme_bw() +
  theme(text = element_text(size = font_size), axis.text = element_text(size = font_size),
        strip.text = element_text(size = 8), legend.position = 'None') + 
  scale_color_discrete(name = "nu") + 
  scale_colour_brewer(palette=7) + 
  #  geom_point(aes(x = cellID, y = rna)) +
  NULL
ggsave("plots/fig2b_without_schematic.eps",device = my_device)
#add schematic image of egg chamber to plot
im <- magick::image_read('plots/fig1c.eps')
df <- data_frame(x = 13,
                 y = 1900,
                 width = 1600,
                 image = list(im))
p3 <- p3 + geom_subview(aes(x=x,y=y,subview=image,width=width,height=width), data=df)
  if (animate_on){
  library(gganimate)
  p3 <- p3 + labs(title = 'Time: {frame_time}', x = "Cell ID", y = "mRNA") +
  transition_time(time) +
  ease_aes('linear')
  }
print(p3)
ggsave('plots/fig2b.eps',device = my_device)

p5 <- ggplot(all_extracted_samples %>%
               filter(nu>0.85 & nu<0.95, cellID==1),
             aes(x = time, y = median, group=factor(nu))) +
  geom_line() +
  labs(title='a)', x = "Time (hrs)", y = "mRNA\ncomplexes") +
  theme_bw() +
  theme(text = element_text(size = font_size), axis.text = element_text(size = font_size),
        strip.text = element_text(size = 8), legend.position = 'None') + 
  scale_color_discrete(name = "nu") + 
  scale_colour_brewer(palette=7) + 
  #  geom_point(aes(x = time, y = rna)) +
  NULL
print(p5)
ggsave('plots/fig2d.eps',device = my_device)

if (!animate_on){
  library(patchwork)
  p3 <- p3 + labs(title='b) Time: 32 hrs')
  p4 <- p2 + p3 + plot_layout(ncol = 1,heights = c(2, 1))
  print(p4)
  ggsave('plots/fig2.eps', width=9,height=9,device = my_device)
  
  library(patchwork)
  p4 <- p5 + p3 + plot_layout(ncol = 1,heights = c(1, 1))
  print(p4)
  ggsave('plots/fig2v2.eps', width=9,height=9,device = my_device)
}

