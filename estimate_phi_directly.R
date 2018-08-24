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
#########################
#process by subtracting background values and valculating ratio
processed <- particles %>% 
  group_by(Sample) %>% 
  mutate(MeanBgd = mean(RawIntDen[Region=='background']),BgdSubtract = RawIntDen-MeanBgd)

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
  out_OE <- estimate_phi_directly(is_wildtype = FALSE, v=seq_len(8)[-3])
  out_WT <- estimate_phi_directly(is_wildtype = TRUE, v=seq_len(11)[-1])
  out_UE <- estimate_phi_directly(is_wildtype = FALSE, is_overexpressor=FALSE, v=seq_len(5))
  a <- out_OE[[1]] %>% mutate(phenotype = 'OE')
  b <- out_WT[[1]] %>% mutate(phenotype = 'WT')
  d <- out_UE[[1]] %>% mutate(phenotype = 'UE')
  h <- full_join(a,b) %>%
    full_join(.,d) %>% ggplot(.,aes(Region,BgdSubtract/MeanByRegion[Region=='nurse_cells'],color=phenotype)) +
     geom_violin(draw_quantiles = c(0.5)) + 
     geom_jitter(alpha=0.3) +
     theme_bw() +
     scale_x_discrete("Region", labels = c("background" = "BGD","oocyte" = "OO","nurse_cells" = "NC")) +
     ylab('Normalised intensity')
  print(h)
  ggsave('plots/pheno_temp.eps',device = cairo_ps)
}

get_mean_and_std <- function(q){
  q %>% mutate(phi = median(q$MeanByRegion) / BgdSubtract) %>%
    summarise(av = mean(phi), std = sd(phi))
}

plot_distn_of_norm_int <- function(){
  get_phi <- function(q) q %>% mutate(phi = median(q$MeanByRegion) / BgdSubtract) %>%
    mutate(phi = 1/phi)
  out_OE <- estimate_phi_directly(is_wildtype = FALSE, v=seq_len(8))
  out_WT <- estimate_phi_directly(is_wildtype = TRUE, v=seq_len(11))
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
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  geom_jitter(aes(color=factor(Stage)),alpha=0.5) + 
  facet_wrap(~phenotype) + 
  coord_cartesian(ylim = c(0, 20)) 
  
}
