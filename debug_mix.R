library(rstan)
expose_stan_functions('debug.stan')
#This constructs a matrix B to illustrate the behaviour. B should be a 16x16 matrix with rank 14
B = t(alter_matrix(construct_matrix(0),c(5,13))) 

Q = qr.Q(qr(B)) #qr decomposition with R
Q_stan = qr_Q_stan(B) #qr decomposition with Stan
apply(Q-Q_stan,2,function(x) sum(abs(x))<10^-14) #several of the final columns are different

#the final two columns of the Q matrix should span the null space of B. If both span the same space, then rank of this should be 2, but gives rank 3.
X = cbind(Q[,15:16],Q_stan[,get_index_for_Q_cols(B)])
print(qr(X)$rank)  
######################################

Q
R = qr.R(qr(B))
diag(R)
R_stan = qr_R_stan(B)
diag(R_stan)

##########################################


expose_stan_functions('model_comparison_at_stst4.stan')

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
rc_indices = seq(from=2,to=16)
names(rc_indices) = seq_len(16)[-1]
p2 <- rc_indices %>%
  purrr::map(get_RC_from_dict) %>%
  purrr::map_df(function(x) get_k2_R(x,nu=nu)) %>%
  mutate(cellID = seq_len(16)) %>%
  gather(key=rcID,value=mRNA,-cellID) %>%
  mutate(rcID=as.integer(rcID)) %>%
  ggplot(aes(x=cellID,y=mRNA)) + 
  geom_line() + 
  facet_wrap(~rcID,scales='free')
print(p1)

my_fun <- function(Q){
  qqr <- Q[,15:16] %>% qr
  qqr$rank
}
