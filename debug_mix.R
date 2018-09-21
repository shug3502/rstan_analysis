library(rstan)
expose_stan_functions('model_comparison_at_stst4.stan')
expose_stan_functions('debug.stan')
nu=0.5
B = alter_matrix(construct_matrix(nu),c(1,5))
Q = qr.Q(qr(t(B)))
Q_stan = qr_Q_stan(t(B))
#Q = Q_stan

compute_N_bar <- function(Q){
N=Q[,15]
N_tilde=Q[,16]

#normalise
N = N/N[min(which(abs(N)>0))]
N_tilde = N_tilde/N_tilde[min(which(abs(N_tilde)>0))]

N_tilde = N - N_tilde;
s_tilde = get_first_nonzero_entry(N_tilde);
N_tilde = my_normalise(N_tilde,s_tilde[1]);
N = N - N_tilde*N[my_floor(s_tilde[2])];


N_bar = (13/15*N + 2/15*sum(N)/sum(N_tilde)*N_tilde)*15/13 
#plot(N_bar)
return(N_bar)
}

compute_N_bar(Q)
compute_N_bar(Q_stan)

################
library(dplyr)
library(tidyr)
library(ggplot2)
nu = 0.95
get_k2_R <- function(x,nu=0.95) {
  B = alter_matrix(construct_matrix(nu),x)
  Q = qr.Q(qr(t(B)))
  return(compute_N_bar(Q))
}
rc_indices = seq_len(16)
names(rc_indices) = seq_len(16)
p1 <- rc_indices %>%
  purrr::map(get_RC_from_dict) %>%
  purrr::map_df(function(x) get_k2_R(x,nu=nu)) %>%
  mutate(cellID = seq_len(16)) %>%
  gather(key=rcID,value=mRNA,-cellID) %>%
  mutate(rcID=as.integer(rcID)) %>%
  ggplot(aes(x=cellID,y=mRNA)) + 
  geom_line() + 
  facet_wrap(~rcID,scales='free')
print(p1)

