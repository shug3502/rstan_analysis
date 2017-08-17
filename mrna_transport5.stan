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
  sigma ~ normal(1.0,0.25) T[0,]; //cauchy(0,2.5) T[0,];
  theta ~ normal(0,10);
  z = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B, x_i);
  for (t in 1:T){
    for (j in 1:16) {
      y[t,j] ~ poisson(sigma*z[t,j]);
    }
  }
}
generated quantities {
  int y_pred[T,16];
  real y_ode[T,16];
  y_ode = integrate_ode_rk45(mrnatransport, y0, t0, ts, theta, B, x_i );
  for (t in 1:T){
    for (j in 1:16){
      y_pred[t,j] = poisson_rng(sigma*y_ode[t,j]);
    }
  }
}
 

 