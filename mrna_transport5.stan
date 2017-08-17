functions {
  real[] mrnatransport(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i
  ) {
    vector[16] dydt;
    matrix[16,16] B;
    B = to_matrix(x_r,16,16);
    dydt = theta[1] * B * to_vector(y) + theta[2];
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
  //  real y0[16];
  real<lower=0> sigma;
  real<lower=0> theta[2];
}
model {
  real z[T,16];
  int cell_indices[5]; //only model these cells
  sigma ~ cauchy(0,2.5) T[0,]; //normal(1.0,0.25) T[0,]; 
  theta ~ normal(0,10);
  z = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B, x_i);
  cell_indices[1] = 1;
  cell_indices[2] = 2;
  cell_indices[3] = 3;
  cell_indices[4] = 9;
  cell_indices[5] = 5;
  for (t in 1:T){
    for (j in 1:5) {
//      y[t,j] ~ poisson(sigma*z[t,j]);
      y[t,j] ~ neg_binomial_2(z[t,cell_indices[j]], sigma);
    }
  }
}
generated quantities {
  int y_pred[T,16];
  real y_ode[T,16];
  y_ode = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B, x_i );
  for (t in 1:T){
    for (j in 1:16){
//      y_pred[t,j] = poisson_rng(sigma*y_ode[t,j]);
      y_pred[t,j] = neg_binomial_2_rng(y_ode[t,j], sigma);
    }
  }
}
 

 