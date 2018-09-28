#mrna transport analytic solution and sensitivity analysis
#could do as markdown
#JH 14/09/17
###########################################

library(Matrix)
library(MASS)
library(limSolve)
library(rstan)
library(ggplot2)
library(dplyr)
source('get_nc_transition_matrix.R')
B1 = get_nc_transition_matrix(0)
B2 = get_nc_transition_matrix(1)
expose_stan_functions(stan_model('model_comparison_at_stst2.stan'))

############################################
#set ICs
#nu = 0.72
th = c(0.24,12.9)
#sig = 100 
#phi = 0.3
#y_0 = rep(0,16)

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

for (ii in c(TRUE,FALSE)){ #make plots for both dyda and dydb

nu = 0.95
t = 30
df <- data.frame(a=numeric(),b=numeric(),dy=numeric())
for (a in seq(from=-3,to=3,length=28)){
  for (b in seq(from=-3,to=3,length=28)){
    temp = c(a=10^a,b=10^b,dy=diff_y(10^a,10^b,nu,t,ii) %>% norm(.,type="1"))
    df <- df %>% bind_rows(.,temp)
  }
}


############################################
g <- ggplot(df, aes(x=log10(a), y=log10(b)))
g <- g + geom_raster(aes(fill=log10(dy)),interpolate=TRUE) +
    geom_contour(aes(z=log10(dy),color='darkorange')) + 
  theme_bw() + 
  theme(text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)) +
  scale_fill_gradientn(colours=c('cyan','dodgerblue4'),limits=c(-8,12))#colours=c('lightpink4','mediumorchid1'))
pts <- data.frame(a=c(log10(th[2])),b=c(log10(th[1])))
g <- g + geom_point(data=pts,aes(a,b),color='darkorange',size=5) +
guides(color='none')
print(g)
plot_name = paste('plots/sensitivity_dyd',ifelse(ii,'a','b'),'.eps',sep='')
ggsave(plot_name,device=cairo_ps)
}

