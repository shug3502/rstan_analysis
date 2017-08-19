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
  int<lower=1> T;
  int y[T,16];
  real y0[16];
  real t0;
  real ts[T];
  real B[16*16];
}
transformed data {
  real x_r[16*16];
  int x_i[0];
}
parameters {
//  real<lower=0> sigma; //noise param
  real<lower=0,upper=1> phi; //difference between particles in NCs and in Oocyte
  real<lower=0> theta[2];
}
model {
  real z[T,16];
  int cell_indices[5]; //only model these cells
  //sigma ~ normal(0,10) T[0,]; //cauchy(0,2.5) T[0,]; //normal(1.0,0.25) T[0,]; 
  phi ~ normal(0,0.1) T[0,1];
  theta ~ cauchy(0,3.5); //normal(0,10);
  z = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B, x_i);
  cell_indices[1] = 1;
  cell_indices[2] = 2;
  cell_indices[3] = 3;
  cell_indices[4] = 9;
  cell_indices[5] = 5;
  for (t in 1:T){
    for (j in 1:5) {
      if (j>1){
        //y[t,j] ~ neg_binomial_2(z[t,cell_indices[j]], sigma);
        y[t,j] ~ poisson(z[t,cell_indices[j]]);
      } else {
        //y[t,j] ~ neg_binomial_2(phi*z[t,cell_indices[j]], sigma);    
        y[t,j] ~ poisson(phi*z[t,cell_indices[j]]);    
      }
    }
  }
}
generated quantities {
  int y_pred[T,16];
  real y_ode[T,16];
  y_ode = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B, x_i );
  for (t in 1:T){
    for (j in 1:16){
      if (j>1){
        //y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j], sigma);
        y_pred[t,j] = poisson_rng(y_ode[t,j]);
      } else {
        //y_pred[t,j] = neg_binomial_2_rng(phi*y_ode[t,j], sigma);    
        y_pred[t,j] = poisson_rng(phi*y_ode[t,j]);    
      }
    }
  }
}
 

 
