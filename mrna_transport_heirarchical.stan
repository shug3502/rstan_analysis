functions {
  real[] mrnatransport(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    vector[16] dydt;
    matrix[16,16] B;
    vector[16] producers;
    producers = rep_vector(1,16);
    producers[1] = 0;
    B = to_matrix(x_r,16,16);
    dydt = theta[1] * B * to_vector(y) + theta[2] * producers;
    return to_array_1d(dydt);
  }
}
data {
  int<lower=1> T1;
  int<lower=1> T2;
  real y[T1,16];
  real y0[16];
  real t0;
  real ts1[T1];
  real ts2[T2];
  real B[16*16];
}
transformed data {
  real x_r[16*16];
  int x_i[0];
}
parameters {
  real<lower=0> sigma; //noise param
  real<lower=0> mu[2];
  //real<lower=0,upper=1> psi;
  real<lower=0> tau;
  //real<lower=0> zeta;
  real<lower=0,upper=1> phi; //difference between particles in NCs and in Oocyte
  real<lower=0> theta[T1,2];
}
model {
  real z[T1,16];
  int cell_indices[5]; //only model these cells
  real aux[T1,16];
  mu ~ normal(0,50);
  tau ~ cauchy(0,2.5);
  //psi ~ normal(0.289,0.0285);
  //zeta ~ normal(0,0.0285);
  phi ~ normal(0.289,0.0285);
  sigma ~ normal(0,50) T[0,]; //cauchy(0,2.5) T[0,]; //normal(1.0,0.25) T[0,]; 
  for (t in 1:T1){
    //phi[t] ~ normal(psi,zeta) T[0,1];
    theta[t] ~ normal(mu,tau);
    aux = integrate_ode_rk45(mrnatransport, y0, t0, ts1, theta[t], B, x_i);
    z[t] = aux[t];
  }
  cell_indices[1] = 1;
  cell_indices[2] = 2;
  cell_indices[3] = 3;
  cell_indices[4] = 9;
  cell_indices[5] = 5;
  for (t in 1:T1){
    for (j in 1:5) {
      if (j>1){
        //y[t,j] ~ neg_binomial_2(z[t,cell_indices[j]], sigma);
        //y[t,j] ~ poisson(z[t,cell_indices[j]]);
      	y[t,j] ~ normal(z[t,cell_indices[j]],sigma);
      } else {
        //y[t,j] ~ neg_binomial_2(phi[t]*z[t,cell_indices[j]], sigma);    
        //y[t,j] ~ poisson(phi[t]*z[t,cell_indices[j]]);    
        y[t,j] ~ normal(phi*z[t,cell_indices[j]],sigma);
      }
    }
  }
}
generated quantities {
  real y_pred[T2,16];
  real y_ode[T2,16];
  real temp[T2,16];
  real theta_pred[T2,2];
  //real phi_pred[T2];
  for (t in 1:T2){ 
    //phi_pred[t] = normal_rng(psi,zeta);
    theta_pred[t,1] = normal_rng(mu[1],tau);
    theta_pred[t,2] = normal_rng(mu[2],tau);
    temp = integrate_ode_rk45(mrnatransport, y0, t0, ts2, theta_pred[t], B, x_i );
    y_ode[t] = temp[t];
    for (j in 1:16){
      if (j>1){
        //y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j], sigma);
        //y_pred[t,j] = poisson_rng(y_ode[t,j]);
      	y_pred[t,j] = normal_rng(y_ode[t,j],sigma);
      } else {
        //y_pred[t,j] = neg_binomial_2_rng(phi[t]*y_ode[t,j], sigma);    
        //y_pred[t,j] = poisson_rng(phi[t]*y_ode[t,j]);    
        y_pred[t,j] = normal_rng(phi*y_ode[t,j],sigma);
      }
    }
  }
}
 

 
