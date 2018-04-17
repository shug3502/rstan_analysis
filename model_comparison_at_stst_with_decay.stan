functions {
  matrix construct_matrix(real th, real gamma_scaled) {
  //add decay at rate gamma
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
  return   B' - diag_matrix(rep_vector(gamma_scaled, 16));
  }
 vector get_k2(real th, real gamma_scaled){
    matrix[16,16] B;
    matrix[16,16] Q;
    vector[16] N;
    real s;
    //int r;
    B = construct_matrix(th, gamma_scaled);    
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
  real<lower=0,upper=1> phi;
  real<lower=0> gamma;
}
model {
  int cell_indices[16]; //only model these cells
  real y_stst[T1,16];
  xi ~ normal(0,0.05) T[0,];
  nu ~ beta(1,1) T[0,1];
  phi ~ normal(0.289,0.0285) T[0,1];
  gamma ~ normal(0,0.01) T[0,];
//for (j in 1:16){ cell_indices[j] = j; }
  cell_indices[1] = 1;
  cell_indices[2] = 2;
  cell_indices[3] = 3;
  cell_indices[4] = 9;
  cell_indices[5] = 5;
for (t in 1:T1){
  //relying on the fact that in practice dim of null space is 1, unless nu=0 (unidirectional backward transport)
  y_stst[t] = to_array_1d(get_k2(nu, gamma));
  for (j in 1:5){
    if (j>1) {
      y_obs[t,cell_indices[j]] ~ normal(y_stst[t,cell_indices[j]]/phi,xi);
    } else {
      y_obs[t,cell_indices[j]] ~ normal(y_stst[t,cell_indices[j]],xi);
    }
  }
}
}
generated quantities {
  real y_pred[(T1+T2),16];
  real y_sim[(T1+T2),16];
  vector[T1] log_lik;
  int cell_indices[16]; //only model these cells
  real y_stst[T1,16];
  cell_indices[1] = 1;
  cell_indices[2] = 2;
  cell_indices[3] = 3;
  cell_indices[4] = 9;
  cell_indices[5] = 5;

  for (t in 1:(T1+T2)) {
    y_pred[t] = to_array_1d(get_k2(nu, gamma));
    for (i in 1:16){
      if (i>1) {
        y_sim[t,i] = normal_rng(y_pred[t,i]/phi,xi);
      } else {
        y_sim[t,i] = normal_rng(y_pred[t,i],xi);
      }
    }
  }

  //compute log likelihood for model comparison via loo
  log_lik = rep_vector(0,T1);
  for (t in 1:T1){
    y_stst[t] = to_array_1d(get_k2(nu, gamma));
    for (j in 1:5){
      if (j>1) {
        log_lik[t] = log_lik[t] + normal_lpdf(y_obs[t,cell_indices[j]] | y_stst[t,cell_indices[j]]/phi,xi);
      } else {
        log_lik[t] = log_lik[t] + normal_lpdf(y_obs[t,cell_indices[j]] | y_stst[t,cell_indices[j]],xi);
      }
    }
  }
}

 
