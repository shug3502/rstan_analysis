library(dplyr)
library(ggplot2)

get_stage_info <- function(root, v){
  stages <- rep(0,length(v)) #don't need to extract estimated stages for each egg chamber example
  for (j in seq_along(v)){
    temp = read.table(file = paste('data/',root,v[j],'/','filenames.txt',sep='')) %>%
      unlist %>%
      stringr::str_extract(., 'stg.') %>%
      stringr::str_split(.,'stg',simplify=TRUE)
    stages[j] = temp[1,2] %>% as.numeric
  }
  return(stages)
}

estimate_phi_directly <- function(is_wildtype, is_overexpressor=TRUE, v=seq(from=2,to=6)){
  #wrap up code to estimate phi based on comparing brightness of spots in different regions of the egg chamber and comparing ratios
  #this allows us to consider whether the same aggregation is occuring in WT and mutant embryos
  #JH 02/07/2018
################################  
  #setup
#v points towards the labelled data
regions <- c('background','nurse_cells','oocyte')
root <- case_when(is_wildtype ~ 'Example',
                  is_overexpressor ~ 'Overexpression',
                  TRUE ~ 'Underexpression')
stages <- get_stage_info(root, v)

particles <- data.frame(Area=numeric(),Mean=numeric(),Min=numeric(),Max=numeric(),RawIntDen=numeric(),
                        Region=factor(),Sample=factor(),Stage=integer())
#####################
#loop over each regions
for (r in seq_along(regions)){
  #get list of files for all examples, and extract results for a given region
  #returns a list of results for all examples for a given region
  temp <- list.files(path = paste('data/',root,v,'/',sep=''), pattern = "\\.xls$",full.names=TRUE) %>%
    grep(regions[r],. , value=TRUE) %>% 
    lapply(.,read.csv) #somehow stored data in different formats
    
  #take results for each example and add to data frame  
  for (j in seq_along(temp)){
    M <- nrow(temp[[j]]) #how many particles are labelled for each region and each egg chamber
    Region <- rep(regions[r],M)
    Sample <- rep(j,M)
    Stage <- rep(stages[j],M)
    particles <- rbind(particles,cbind(cbind(cbind(temp[[j]],Region),Sample),Stage))  
  }
}
gg<- particles %>% filter(Region=="background") %>% ggplot(aes(x=RawIntDen, y=Region, color=factor(Sample))) + geom_jitter()
print(gg)
#########################
#process by subtracting background values and valculating ratio
processed <- particles %>% 
  filter(Region!="background" | RawIntDen<10000) %>% 
  group_by(Sample) %>% 
  mutate(MeanBgd = median(RawIntDen[Region=='background']),BgdSubtract = RawIntDen-MeanBgd)

q <- processed %>%
  dplyr::select(BgdSubtract,Region,Sample,Stage) %>% #be careful of clash with select in MASS package
  group_by(Region) %>%
  mutate(MeanByRegion = mean(BgdSubtract),StdByRegion = sd(BgdSubtract))
q <- q %>%
  mutate(Ratio = MeanByRegion/median(q$MeanByRegion))
hh <- ggplot(q,aes(Region,BgdSubtract/MeanByRegion[Region=='nurse_cells'],color=factor(Sample))) +
  geom_violin() + 
  geom_jitter() +
  theme_bw() +
  scale_x_discrete("Region", labels = c("background" = "BGD","oocyte" = "OO","nurse_cells" = "NC")) +
  ylab('Normalised intensity')

return(list(q,hh))
}

make_plot_comparing_phenotypes <- function(){
  out_OE <- estimate_phi_directly(is_wildtype = FALSE, v=seq_len(13))
  out_WT <- estimate_phi_directly(is_wildtype = TRUE, v=seq_len(12))
  out_UE <- estimate_phi_directly(is_wildtype = FALSE, is_overexpressor=FALSE, v=seq_len(5)[-1])
  a <- out_OE[[1]] %>% mutate(phenotype = 'OE')
  b <- out_WT[[1]] %>% mutate(phenotype = 'WT')
  d <- out_UE[[1]] %>% mutate(phenotype = 'UE')
  h <- full_join(a,b) %>%
    full_join(.,d) %>% ggplot(.,aes(Region,BgdSubtract/MeanByRegion[Region=='nurse_cells'],color=phenotype)) +
     geom_violin(draw_quantiles = c(0.5), scale="width") + 
     geom_jitter(alpha=0.1) +
     theme_bw() +
     scale_x_discrete("Region", labels = c("background" = "BGD","oocyte" = "OO","nurse_cells" = "NC")) +
     ylab('Normalised intensity')
  print(h)
  ggsave('plots/pheno_temp.eps',device = cairo_ps)
}

process_for_phi <- function(q){
  z<- q %>% 
    mutate(phi = median(q$MeanByRegion) / BgdSubtract) %>%
    ungroup() %>%
    filter(Region=="oocyte") %>%
    filter(phi>-1 & phi<10) %>%
    group_by(phenotype)
  return(z)
}

get_mean_and_std <- function(q){
 q %>%
    process_for_phi() %>%
    add_tally() %>%
    summarise(av = mean(phi), std = sd(phi), av_median=median(phi), stnd_error=sd(phi)/n[1])
}

#attempt to fit gamma distribution to use as prior for phi
fit_gamma_to_phi_data <- function(q){
  library(fitdistrplus)
  z <- q %>%
    process_for_phi() %>% 
    summarise(shape = fitdist(phi, distr = "gamma", method = "mle")$estimate[1], rate = fitdist(phi, distr = "gamma", method = "mle")$estimate[2])
  # fit.gamma <- fitdist(z$phi, distr = "gamma", method = "mle")
  # summary(fit.gamma)
  # plot(fit.gamma)
  
  return(z)
}


plot_distn_of_norm_int <- function(){
  get_phi <- function(q) q %>% mutate(phi = median(q$MeanByRegion) / BgdSubtract) %>%
    mutate(phi = 1/phi)
  out_OE <- estimate_phi_directly(is_wildtype = FALSE, v=seq_len(13))
  out_WT <- estimate_phi_directly(is_wildtype = TRUE, v=seq_len(12))
  out_UE <- estimate_phi_directly(is_wildtype = FALSE, is_overexpressor=FALSE, v=seq_len(5)[-1])
  a <- out_OE[[1]] %>% mutate(phenotype = 'OE')
  b <- out_WT[[1]] %>% mutate(phenotype = 'WT')
  d <- out_UE[[1]] %>% mutate(phenotype = 'UE')
  
  full_join(a,b) %>%
  full_join(.,d) %>%
  get_phi %>%
  filter(Region=='oocyte') %>% 
#      filter(phi>0 ) %>%
  ggplot(aes(x = Region,y = phi)) + 
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), scale="width") +
  geom_jitter(aes(color=factor(Sample)),alpha=0.5) + 
  facet_wrap(~phenotype) + 
  coord_cartesian(ylim = c(0, 20)) 
  
}

out_OE <- estimate_phi_directly(is_wildtype = FALSE, v=seq_len(13))
out_WT <- estimate_phi_directly(is_wildtype = TRUE, v=seq_len(12))
out_UE <- estimate_phi_directly(is_wildtype = FALSE, is_overexpressor=FALSE, v=seq_len(5)[-1])
a <- out_OE[[1]] %>% mutate(phenotype = 'OE')
b <- out_WT[[1]] %>% mutate(phenotype = 'WT')
d <- out_UE[[1]] %>% mutate(phenotype = 'UE')
q <- full_join(a,b) %>%
  full_join(.,d)
get_mean_and_std(q) %>% print()
fit_gamma_to_phi_data(q) %>% print()

q %>%
  process_for_phi() %>%
  group_by(phenotype) %>% 
  ggplot(aes(phi,color=phenotype)) +
  geom_density() 

q %>% process_for_phi() %>%
  group_by(Sample,phenotype) %>%
  summarise(av=mean(phi)) %>%
  ggplot(aes(x=phenotype,y=av)) +
  geom_violin(draw_quantiles = c(0.5), scale="width") +
  geom_jitter() + 
  labs(y="phi") + 
  theme_bw()

q %>% process_for_phi() %>%
  group_by(Sample,phenotype) %>%
  summarise(av=mean(phi)) %>%
  group_by(phenotype) %>%
  add_tally() %>%
  summarise(ss=mean(av),std_err=sd(av)/sqrt(n[1]))

a %>% 
  mutate(Normalised_intensity = BgdSubtract / median(a$MeanByRegion)) %>%
  ggplot(aes(x=Region,y=Normalised_intensity)) +
  geom_violin(draw_quantiles = c(0.5)) +
  geom_jitter(alpha=0.3) +
  labs(y='Normalised intensity') +
  scale_x_discrete("Region", labels = c("background" = "BGD","oocyte" = "OO","nurse_cells" = "NC")) +
  theme_bw()
ggsave('plots/fig3d.eps',device=cairo_ps)

a %>% 
  mutate(Normalised_intensity = BgdSubtract / median(a$MeanByRegion)) %>%
  filter(Region=='oocyte') %>%
  ggplot(aes(Normalised_intensity)) +
  geom_density() +
  labs(x='Normalised intensity', y='Density') +
#  scale_x_discrete("Region", labels = c("background" = "BGD","oocyte" = "OO","nurse_cells" = "NC")) +
  theme_bw()

#distribution is not normal as shown by a density plot
a %>% 
  mutate(Normalised_intensity = BgdSubtract / median(a$MeanByRegion)) %>%
  summarise(ava=mean(Normalised_intensity),meda=median(Normalised_intensity))

a %>% 
  mutate(phi = median(a$MeanByRegion)/BgdSubtract) %>%
  summarise(ava=mean(phi),meda=median(phi))

