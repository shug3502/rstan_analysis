library(rstan)
library(dplyr)
expose_stan_functions('model_comparison_at_stst4.stan') #for functions to get k2 and B etc

analytic_soln <- function(b,a,nu,t,y0=rep(0,16),producers=c(0,rep(1,15))){
  #implement analytic solution as a function of model parameters
  #b,a,nu,t should be numeric scalar values
  ############
  B = construct_matrix(nu)
  k2 = get_k2(nu,0)
  k2 = k2/sum(k2) #normalise so that the sum gives 1
  k1 = solve(-b*B[2:16,2:16],a*producers[2:16]-k2[2:16])
  k1 = c(1,k1)
  e = eigen(b*B)
  D = diag(exp(e$values*t))
  V = e$vectors
  x = solve(V,y0-k1)
  y = V%*%D%*%x + k1 + 15*a*k2*t
  return(y)
}


observation_model <- function(phi,sigma,y){
  z = y
  z[1] = y[1]*phi
  #put negative binomial stuff in here with correct parameterisation
  return(z)
}
