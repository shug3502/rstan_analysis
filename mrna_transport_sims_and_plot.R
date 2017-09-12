setwd('~/Documents/FISH_data/rstan_analysis')
library(rstan)
library(mvtnorm)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
identifier = 'v002' #run identifier
nSamples = 15 #how many egg chambers segmented

#############################################################
#Fit linear model to log(length) of egg chambers to get time of development.
#Could use area or volume of egg chamber. Compare to Jia 2016 Sci Reports.
#Area may be better as in some examples egg chambers are quite squashed.
#To get lengths read length.txt file in each folder

egg_chamber_areas <- rep(0,nSamples)
#stages <- rep(0,nSamples) #don't need to extract estimated stages for each egg chamber example
for (j in 1:nSamples){
  egg_chamber_areas[j] <- as.numeric(read.table(paste('../data/Example',j,'/area.txt',sep='')))
  #   temp <- list.files(path = paste('../data/Example',j,'/',sep=''), pattern = "\\_grk\\.tif$") %>%
  #     stringr::str_extract(., 'stg.') %>%
  #     stringr::str_split(.,'stg',simplify=TRUE)
  #   stages[j] = temp[2] %>% as.numeric
}
#lm_time <- lm(log(egg_chamber_lengths) ~ stages)

#take l0 = log(20) as initial time (when using length)
#take measure egg chamber lengths as a scaled time variable
t0 = log(400)
log_areas = sort.int(log(egg_chamber_areas),index.return=TRUE)
ts = log_areas$x #ts needs to be an ordered time vector
sort_indices = log_areas$ix

#############################################################
m0 = c(0, rep(1,15)) #initial condition
th = c(2,2)
sig = 10^-9
phi = 1
deltaT = 0.1
t0 = 0.0
ts = seq(deltaT,nSamples * deltaT,deltaT)

#############################################################
source('get_nc_transition_matrix.R')
B1 = get_nc_transition_matrix(0) %>% as.vector
B2 = get_nc_transition_matrix(1) %>% as.vector

#############################################################

  raw_data = matrix(c(1840.,   283.,   194.,   168.,   190.,   179.,    67.,    22.,
     85.,    70.,    56.,    47.,   174.,   136.,    20.,    17.,
     708.,    47.,   117.,    50.,    91.,    34.,    86.,    49.,
      76.,    94.,   100.,    63.,   110.,    73.,    81.,    71.,
     1454.,   261.,   144.,   180.,   108.,   142.,    10.,    96.,
     18.,   127.,    20.,    53.,    19.,   109.,    20.,    11.,
     1106.,   162.,    36.,   118.,    92.,   131.,    61.,    76.,
     82.,   129.,    39.,   135.,   122.,    84.,    31.,    35.,
    1485.,   407.,    53.,    NA,   399.,    NA,    NA,    NA,
       21.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
      1652.,    68.,   376.,    NA,    43.,    NA,    NA,    NA,
       263.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
      1012.,    23.,   141.,    NA,    28.,    NA,    NA,    NA,
       88.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
       800.,    44.,   119.,    NA,    29.,    NA,    NA,    NA,
        5.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
      1007.,    51.,     5.,    NA,    32.,    NA,    NA,    NA,
       11.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
      1045.,    78.,   130.,    NA,   101.,    NA,    NA,    NA,
       11.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
       446.,    83.,    22.,    NA,     2.,    NA,    NA,    NA,
        44.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
     1862.,    48.,    32.,    NA,    73.,    NA,    NA,    NA,
      21.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
      566.,    13.,    18.,    NA,    16.,    NA,    NA,    NA,
       10.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
     1032.,    31.,   144.,    NA,    43.,    NA,    NA,    NA,
      6.,    NA,    NA,    NA,    NA,    NA,    NA,    NA,
     1076.,    22.,    19.,    NA,    67.,    NA,    NA,    NA,
      70.,    NA,    NA,    NA,    NA,    NA,    NA,    NA ),ncol=16,byrow=TRUE)
  raw_data = raw_data[sort_indices,] #need to sort time series and correspondingly reorder rows

  #sample from the model to get fake data
  print('using fake generated data')
  samples <- stan(file = 'mrna_transport6.stan',
                  data = list (
                    T  = nSamples,
                    y0 = m0,
                    t0 = t0,
                    ts = ts,
                    theta = array(th, dim = 2),
                    sigma = sig,
                    phi = phi,
                    B = B2
                  ),
                  algorithm="Fixed_param",
                  seed = 42,
                  chains = 1,
                  iter =100, 
                  refresh = -1
  )
  
 # exp_data = matrix(s[1,1,1:(16*nSamples)],nrow=nSamples,byrow=FALSE) #this is our fake data

source('post_pred_plot.R')
post_pred_plot(raw_data=NA,ts,nSamples,'y_hat',samples,identifier,title_stem='plots/sims')

