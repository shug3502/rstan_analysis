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
    //int r;
    B = construct_matrix(th);    
    //r = 15; //compute_rank(B); //rank will be 15 unless nu=0
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
  int<lower=1> T2;
  real y_obs[T1,16];
}
/*
transformed data {
  real y[T1,16];
  for (t in 1:T1){
    for (i in 1:16){
      y[t,i] = y_obs[t,i]/y_obs[t,1]; //normalise the vectors for comparison
    }
  }
}
*/
parameters {
  real<lower=0,upper=1> nu[T1];
  real<lower=0,upper=1> mu;
  real<lower=0> tau;
  real<lower=0> sigma;
  //real<lower=0> zeta;
  //real<lower=0> xi;
}
model {
  real y_stst[T1,16];
  mu ~ beta(0.5,0.5) T[0,1];
  tau ~ normal(0,0.1) T[0,];
  //zeta ~ normal(0,0.1) T[0,];
  //xi ~ normal(0,0.1) T[0,];
  sigma ~ normal(0,0.05) T[0,];
for (t in 1:T1){
  nu[t] ~ normal(mu,tau) T[0,1];
  //sigma[t] ~ normal(zeta,xi) T[0,];
//relying on the fact that in practice dim of null space is 1, unless nu=0 (unidirectional backward transport)
  y_stst[t] = to_array_1d(get_k2(nu[t]));
/*
for (i in 1:16) {
    y_stst[t,i] = y_stst[t,i]/y_stst[t,1]; //normalise for comparison
  }
*/
  y_obs[t] ~ normal(y_stst[t],sigma);
}
}
generated quantities {
  real y_pred[T2,16];
  real y_sim[T2,16];
  real nu_sim[T2];
  //real sigma_sim[T2];
for (t in 1:T2) {
  nu_sim[t] = normal_rng(mu,tau); 
  //sigma_sim[t] = fabs(normal_rng(zeta,xi)); //scale parameter should be +ve
  y_pred[t] = to_array_1d(get_k2(nu_sim[t]));
for (i in 1:16){
//    y_pred[t,i] = y_pred[t,i]/y_pred[t,1];
    y_sim[t,i] = normal_rng(y_pred[t,i],sigma);
  }
}
}

 
