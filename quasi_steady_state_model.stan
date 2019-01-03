functions {
  matrix construct_matrix(real th) {
  real nu;
  matrix[16,16] B;
  nu = (1+th)/2;
  B = rep_matrix(0,16,16);
  B[1,1] = -4*(1-nu);
  B[2,1] = nu;
  B[3,1] = nu;
  B[5,1] = nu;
  B[9,1] = nu;
  B[1,2] = 1 - nu;
  B[2,2] = -3*(1-nu) - nu;
  B[4,2] = nu;
  B[6,2] = nu;
  B[10,2] = nu;
  B[1,3] = 1-nu;
  B[3,3] = -2*(1-nu) - nu;
  B[7,3] = nu;
  B[11,3] = nu;
  B[2,4] = 1-nu;
  B[4,4] = -2*(1-nu) - nu;
  B[8,4] = nu;
  B[12,4] = nu;
  B[1,5] = 1-nu;
  B[5,5] = -(1-nu) - nu;
  B[13,5] = nu;
  B[2,6] = 1-nu;
  B[6,6] = -(1-nu) - nu;
  B[14,6] = nu;
  B[3,7] = 1-nu;
  B[7,7] = -(1-nu) - nu;
  B[15,7] = nu;
  B[4,8] = 1-nu;
  B[8,8] = -(1-nu) - nu;
  B[16,8] = nu;
  B[1,9] = 1-nu;
  B[9,9] = -nu;
  B[2,10] = 1-nu;
  B[10,10] = -nu;
  B[3,11] = 1-nu;
  B[11,11] = -nu;
  B[4,12] = 1-nu;
  B[12,12] = -nu;
  B[5,13] = 1-nu;
  B[13,13] = -nu;
  B[6,14] = 1-nu;
  B[14,14] = -nu;
  B[7,15] = 1-nu;
  B[15,15] = -nu;
  B[8,16] = 1-nu;
  B[16,16] = -nu;
//need to return transpose
  return B';
  }

 vector get_k2(real th){
    matrix[16,16] B;
    matrix[16,16] Q;
    vector[16] N;
    real s;
    B = construct_matrix(th); 
    Q = qr_Q(B'); //compute qr decomposition
  //take last n-r cols of Q as basis of null space
  //may need to also exclude case where null space has dim 0, ie r=16
  N = Q[1:16,16];
  s = N[1];
  for (i in 1:16){ //normalise
      N[i] = N[i]/s;
  }
  return N;
  }
}
data {
  int<lower=1> T1;
  int<lower=0> T2;
  real y_obs[T1,16];
}
parameters {
  real<lower=0,upper=1> nu;
  real<lower=0> xi;
  real<lower=0> phi;
}
model {
  real y_stst[T1,16];
  xi ~ normal(0,0.1) T[0,];
  nu ~ beta(1,1) T[0,1];
  phi ~ normal(0.345,0.047); 
for (t in 1:T1){
  //relying on the fact that in practice dim of null space is 1, unless nu=0 (unidirectional backward transport)
  y_stst[t] = to_array_1d(get_k2(nu));
  for (j in 2:16){
    y_obs[t,j] ~ normal(y_stst[t,j]/phi,xi) T[0,];
  }
}
}
generated quantities {
  real y_pred[(T1+T2),16];
  real y_sim[(T1+T2),16];
  vector[T1] log_lik;
  real y_stst_pred[T1,16];
  for (t in 1:(T1+T2)) {
    y_pred[t] = to_array_1d(get_k2(nu));
    y_sim[t,1] = 1; 
    for (j in 2:16){
      y_sim[t,j] = fabs(normal_rng(y_pred[t,j]/phi, xi));
    }
  }
    //compute log likelihood for model comparison via loo
  log_lik = rep_vector(0,T1);
  for (t in 1:T1){
    y_stst_pred[t] = to_array_1d(get_k2(nu));
    for (j in 2:16){
      log_lik[t] = log_lik[t] + normal_lpdf(y_obs[t,j] | y_stst_pred[t,j]/phi, xi)/(1-normal_lcdf(0 | y_stst_pred[t,j]/phi, xi)); //denominator due to truncation
    }
  }
}
