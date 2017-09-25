#mrna transport analytic solution and sensitivity analysis
#could do as markdown
#JH 14/09/17
###########################################

library(Matrix)
library(MASS)
library(limSolve)
library(rstan)
library(dplyr)
source('get_nc_transition_matrix.R')
B1 = get_nc_transition_matrix(0)
B2 = get_nc_transition_matrix(1)
expose_stan_functions(stan_model('model_comparison_at_stst.stan'))

############################################
#set ICs
#nu = 0.72
th = c(5,100)
sig = 100 
phi = 0.289
y_0 = rep(0,16)

############################################
#try to construct as functions 

k1 <- function(a,b,nu){
  #NB RHS is linear in a, LHS linear in b, use this to give derivatives
  v = c(0,rep(1,15))
  k2 = get_k2(nu)/sum(get_k2(nu))*a*sum(v) #normalised
  k = Solve(-b*construct_matrix(nu), a*v - k2) # uses Solve with pseudo inverse since system is underdetermined
  return(k)
}

diff_y <- function(a,b,nu,t,diff_wrt_a = TRUE){
#calculate derivative of y wrt parameters a, b based on analytic results
  B = construct_matrix(nu) #rate matrix
  X = eigen(B) #compute eigenvalues
  V = X$vectors
  D = diag(exp(t*X$values))
  v = c(0,rep(1,15))
  k2 = get_k2(nu)/sum(get_k2(nu))*a*sum(v) #normalised

  dk1dab <- ifelse(rep(diff_wrt_a,16),(1/a)*k1(a,b,nu),-(1/b)*k1(a,b,nu))
  result <- (diag(16) - V%*%D%*%ginv(V))%*%dk1dab + ifelse(diff_wrt_a,k2*t,0)
  return(result)
}
############################################

nu = 0.72
t = 10
df <- data.frame(a=numeric(),b=numeric(),dy=numeric())
for (a in seq(from=-3,to=3,length=28)){
  for (b in seq(from=-3,to=3,length=28)){
    temp = c(a=10^a,b=10^b,dy=diff_y(10^a,10^b,nu,t,FALSE) %>% norm(.,type="1"))
    df <- df %>% bind_rows(.,temp)
  }
}


############################################
library(ggplot2)
g <- ggplot(df, aes(log10(a), log10(b)))
g + geom_raster(aes(fill=log10(dy)),interpolate=FALSE) +
  theme_bw() + 
  scale_fill_gradientn(colours=c('darkorange','dodgerblue2')) #c('lightpink4','mediumorchid1'))




