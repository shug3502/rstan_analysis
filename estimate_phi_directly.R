library(dplyr)
library(ggplot2)

estimate_phi_directly <- function(is_wildtype, is_overexpressor=TRUE, v=seq(from=2,to=6)){
  #wrap up code to estimate phi based on comparing brightness of spots in different regions of the egg chamber and comparing ratios
  #this allows us to consider whether the same aggregation is occuring in WT and mutant embryos
  #JH 02/07/2018
################################  
  #setup
#v <- seq(from=2,to=6) #v points towards the labelled data
regions <- c('background','nurse_cells','oocyte')
particles <- data.frame(Area=numeric(),Mean=numeric(),Min=numeric(),Max=numeric(),RawIntDen=numeric(),Region=factor(),Sample=factor())

#####################
#loop over each regions
for (r in seq_along(regions)){
  #get list of files for all examples, and extract results for a given region
  #returns a list of results for all examples for a given region
  if (is_wildtype){
    root <- 'Example'
    temp <- list.files(path = paste('data/',root,v,'/',sep=''), pattern = "\\.xls$",full.names=TRUE) %>%
      grep(regions[r],. , value=TRUE) %>% 
      lapply(.,read.table)
  } else {
  root <- ifelse(is_overexpressor,'Overexpression','Underexpression')
  temp <- list.files(path = paste('data/',root,v,'/',sep=''), pattern = "\\.xls$",full.names=TRUE) %>%
    grep(regions[r],. , value=TRUE) %>% 
    lapply(.,read.csv) #somehow stored data in different formats
  }  
  #take results for each example and add to data frame  
  for (j in seq_along(temp)){
    Region <- rep(regions[r],nrow(temp[[j]]))
    Sample <- rep(j,nrow(temp[[j]]))
    particles <- rbind(particles,cbind(cbind(temp[[j]],Region),Sample))  
  }
}
#########################
#process by subtracting background values and valculating ratio
processed <- particles %>% 
  group_by(Sample) %>% 
  mutate(MeanBgd = mean(RawIntDen[Region=='background']),BgdSubtract = RawIntDen-MeanBgd)

q <- processed %>%
  dplyr::select(BgdSubtract,Region,Sample) %>% #be careful of clash with select in MASS package
  group_by(Region) %>%
  mutate(MeanByRegion = mean(BgdSubtract),StdByRegion = sd(BgdSubtract))
q <- q %>%
  mutate(Ratio = MeanByRegion/median(q$MeanByRegion))
hh <- ggplot(q,aes(Region,BgdSubtract/MeanByRegion[Region=='nurse_cells'])) +
  geom_violin() + 
  geom_jitter() +
  theme_bw() +
  scale_x_discrete("Region", labels = c("background" = "BGD","oocyte" = "OO","nurse_cells" = "NC")) +
  ylab('Normalised intensity')

return(list(q,hh))
}

make_plot_comparing_phenotypes <- function(){
  out_OE <- estimate_phi_directly(is_wildtype = FALSE, v=c(1,2,5))
  out_WT <- estimate_phi_directly(is_wildtype = TRUE, v=seq(from=2,to=6))
  out_UE <- estimate_phi_directly(is_wildtype = FALSE, is_overexpressor=FALSE, v=c(1))
  a <- out_OE[[1]] %>% mutate(phenotype = 'OE')
  b <- out_WT[[1]] %>% mutate(phenotype = 'WT')
  d <- out_UE[[1]] %>% mutate(phenotype = 'UE')
  full_join(a,b) %>%
    full_join(.,d) %>% ggplot(.,aes(Region,BgdSubtract/MeanByRegion[Region=='nurse_cells'],color=phenotype)) +
     geom_violin() + 
     geom_jitter() +
     theme_bw() +
     scale_x_discrete("Region", labels = c("background" = "BGD","oocyte" = "OO","nurse_cells" = "NC")) +
     ylab('Normalised intensity')
  ggsave('pheno_temp.eps')
}
